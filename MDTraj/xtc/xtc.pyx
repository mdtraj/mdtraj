# cython: c_string_type=str, c_string_encoding=ascii
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################


###############################################################################
# Imports
###############################################################################

import os
import warnings
import cython
cimport cython
import numpy as np
cimport numpy as np
np.import_array()
from mdtraj.utils.arrays import ensure_type
cimport xdrlib

###############################################################################
# globals
###############################################################################

cdef int _EXDROK = 0             # OK
cdef int _EXDRENDOFFILE = 11     # End of file
_EXDR_ERROR_MESSAGES = {
    1: "Header",
    2: "String",
    3: "Double",
    4: "Integer",
    5: "Float",
    6: "Unsigned integer",
    7: "Compressed 3d coordinate",
    8: "Closing file",
    9: " Magic number",
    10: 'Not enough memory',
    12: "File not found"
}

# numpy variable types include the specific numpy of bytes of each, but the c
# variables in our interface file don't. this could get bad if we're on a wierd
# machine, so lets make sure first
if sizeof(int) != sizeof(np.int32_t):
    raise RuntimeError('Integers on your compiler are not 32 bits. This is not good.')
if sizeof(float) != sizeof(np.float32_t):
    raise RuntimeError('Floats on your compiler are not 32 bits. This is not good')

###############################################################################
# Classes
###############################################################################

cdef class XTCTrajectoryFile:
    """XTCTrajectoryFile(filenamee, mode='r', force_overwrite=True, **kwargs)

    Interface for reading and writing to a GROMACS XTC file.
    This is a file-like objec that supports both reading and writing.
    It also supports the context manager ptorocol, so you can use it
    with the python 'with' statement.

    The conventional units in the XTC file are nanometers and picoseconds.
    The format only supports saving coordinates, the time, the md step,
    and the unit cell parametrs (box vectors)

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already exists on disk, should we overwrite it?

    Other Parameters
    ----------------
    min_chunk_size : int, default=100
        In read mode, we need to allocate a buffer in which to store the data
        without knowing how many frames are in the file. This parameter is the
        minimum size of the buffer to allocate.
    chunk_size_multiplier : int, default=1.5
        In read mode, we need to allocate a buffer in which to store the data without knowing how many frames are in
        the file. We can *guess* this information based on the size of the file on disk, but it's not perfect. This
        parameter inflates the guess by a multiplicative factor.

    Examples
    --------
    >>> # read the data from from an XTC file
    >>> with XTCTrajectoryFile('traj.xtc') as f:
    >>>    xyz, time, step, box = f.read()

    >>> # write some random coordinates to an XTC file
    >>> with XTCTrajectoryFile('output.xtc', 'w') as f:
    >>>     f.write(np.random.randn(10,1,3))

    See Also
    --------
    mdtraj.load_xtc : High-level wrapper that returns a ``md.Trajectory``
    """
    cdef xdrlib.XDRFILE* fh
    cdef int n_atoms          # number of atoms in the file
    cdef int frame_counter    # current position in the file, in read mode
    cdef int is_open          # is the file handle currently open?
    cdef int approx_n_frames  # appriximate number of frames in the file, as guessed based on its size
    cdef char* mode           # mode in which the file is open, either 'r' or 'w'
    cdef int min_chunk_size
    cdef float chunk_size_multiplier
    cdef readonly char* distance_unit


    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True, **kwargs):
        """Open a GROMACS XTC file for reading/writing.
        """
        self.distance_unit = 'nanometers'
        self.is_open = False
        self.frame_counter = 0

        if str(mode) == 'r':
            self.n_atoms = 0
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            xdrlib.read_xtc_natoms(filename, &self.n_atoms)
            if self.n_atoms <= 0:
                raise IOError('Malformed XTC file. Number of atoms <= 0. '
                              'Are you sure this is a valid GROMACS xtc file?')

            self.fh = xdrlib.xdrfile_open(filename, b'r')
            if self.fh is NULL:
                raise IOError('File not found: "%s"' % filename)
            self.approx_n_frames = self._estimate_n_frames_from_filesize(os.stat(filename).st_size)

            self.min_chunk_size = max(kwargs.pop('min_chunk_size', 100), 1)
            self.chunk_size_multiplier = max(kwargs.pop('chunk_size_multiplier', 1.5), 0.01)


        elif str(mode) == 'w':
            if force_overwrite and os.path.exists(filename):
                os.unlink(filename)

            self.fh = xdrlib.xdrfile_open(filename, 'w')
            if self.fh is NULL:
                raise IOError('Unable to open file "%s"' % filename)
        else:
            raise ValueError('mode must be one of "r" or "w". '
                             'you supplied %s' % mode)

        for key in kwargs.keys():
            warnings.warn('kwarg "%s" was not recognized or processed' % key)

        self.is_open = True
        self.mode = mode

    def _estimate_n_frames_from_filesize(self, filesize):
        # model: size(bytes) = coefs_[0] * n_frames + coefs_[1]*n_atoms
        #                       + coefs_[2] * n_frames * n_atoms
        #                       + intercept
        # fit on a small training set with a few hundred frames
        coefs_ = [9.93733050e+01,  -6.49891780e-02,   4.74462831e+00]
        intercept_ = 5

        approx_n_frames = (filesize - intercept_ -
                           coefs_[1]*self.n_atoms) / (coefs_[2] * self.n_atoms +
                                                      coefs_[0])

        return approx_n_frames


    def __dealloc__(self):
        self.close()


    def close(self):
        "Close the XTC file handle"
        if self.is_open:
            xdrlib.xdrfile_close(self.fh)
            self.is_open = False

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """read(n_frames=None, stride=None, atom_indices=None)

        Read data from an XTC file

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : int, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it required
            an extra copy, but will save memory.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
            The cartesian coordinates, in nanometers
        time : np.ndarray, shape=(n_frames), dtype=np.float32
            The simulation time, in picoseconds, corresponding to each frame
        step : np.ndarray, shape=(n_frames), dtype=np.int32
            The step in the simulation corresponding to each frame
        box : np.ndarray, shape=(n_frames, 3, 3), dtype=np.float32
            The box vectors in each frame.
        """
        if not str(self.mode) == 'r':
            raise ValueError('read() is only available when file is opened in mode="r"')
        if not self.is_open:
            raise IOError('file must be open to read from it.')

        if n_frames is not None:
            # if they supply the number of frames they want, that's easy
            if not int(n_frames) == n_frames:
                raise ValueError('n_frames must be an int, you supplied "%s"' % n_frames)
            xyz, time, step, box = self._read(int(n_frames), atom_indices)
            return xyz[::stride], time[::stride], step[::stride], box[::stride]


        # if they want ALL of the remaining frames, we need to guess at the chunk
        # size, and then check the exit status to make sure we're really at the EOF
        all_xyz, all_time, all_step, all_box = [], [], [], []

        while True:
            # guess the size of the chunk to read, based on how many frames we think are in the file
            # and how many we've currently read
            chunk = max(abs(int((self.approx_n_frames - self.frame_counter) * self.chunk_size_multiplier)),
                        self.min_chunk_size)

            xyz, time, step, box = self._read(chunk, atom_indices)
            if len(xyz) <= 0:
                break

            all_xyz.append(xyz)
            all_time.append(time)
            all_step.append(step)
            all_box.append(box)

        if len(all_xyz) == 0:
            return np.array([]), np.array([]), np.array([]), np.array([])
        return (np.concatenate(all_xyz)[::stride],
                np.concatenate(all_time)[::stride],
                np.concatenate(all_step)[::stride],
                np.concatenate(all_box)[::stride])


    def _read(self, int n_frames, atom_indices):
        """Read a specified number of XTC frames from the buffer"""

        cdef int i = 0
        cdef int status = _EXDROK
        cdef int n_atoms_to_read

        if atom_indices is None:
            n_atoms_to_read = self.n_atoms
        elif isinstance(atom_indices, slice):
            n_atoms_to_read = len(np.arange(self.n_atoms)[atom_indices])
        else:
            atom_indices = np.asarray(atom_indices)
            if min(atom_indices) < 0:
                raise ValueError('atom_indices should be zero indexed. you gave an index less than zero')
            if max(atom_indices) >= self.n_atoms:
                raise ValueError('atom indices should be zero indexed. you gave an index bigger than the number of atoms')
            n_atoms_to_read = len(atom_indices)

        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz = \
            np.empty((n_frames, n_atoms_to_read, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time = \
            np.empty((n_frames), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.int32_t, mode='c'] step = \
            np.empty((n_frames), dtype=np.int32)
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box = \
            np.empty((n_frames, 3, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] prec = \
            np.empty((n_frames), dtype=np.float32)

        # only used if atom_indices is given
        cdef np.ndarray[dtype=np.float32_t, ndim=2] framebuffer = np.zeros((self.n_atoms, 3), dtype=np.float32)


        while (i < n_frames) and (status != _EXDRENDOFFILE):
            if atom_indices is None:
                status = xdrlib.read_xtc(self.fh, self.n_atoms, <int*> &step[i],
                                         &time[i], <xdrlib.matrix>&box[i,0,0], <xdrlib.rvec*>&xyz[i,0,0], &prec[i])
            else:
                status = xdrlib.read_xtc(self.fh, self.n_atoms, <int*> &step[i],
                                         &time[i], <xdrlib.matrix>&box[i,0,0], <xdrlib.rvec*>&framebuffer[0,0], &prec[i])
                xyz[i, :, :] = framebuffer[atom_indices, :]

            if status != _EXDRENDOFFILE and status != _EXDROK:
                raise RuntimeError('XTC read error: %s' % _EXDR_ERROR_MESSAGES.get(status, 'unknown'))
            i += 1

        if status == _EXDRENDOFFILE:
            xyz = xyz[:i-1]
            box = box[:i-1]
            time = time[:i-1]
            step = step[:i-1]

        self.frame_counter += i

        return xyz, time, step, box


    def write(self, xyz, time=None, step=None, box=None):
        """write(xyz, time=None, step=None, box=None)

        Write data to an XTC file

        Parameters
        ----------
        xyz : np.ndarray, dtype=np.float32, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms, in nanometers
        time : np.ndarray, dtype=float32, shape=(n_frames), optional
            The simulation time corresponding to each frame, in picoseconds.
            If not supplied, the numbers 0..n_frames will be written.
        step :  np.ndarray, dtype=int32, shape=(n_frames), optional
            The simulation timestep corresponding to each frame, in steps.
            If not supplied, the numbers 0..n_frames will be written
        box : np.ndarray, dtype=float32, shape=(n_frames, 3, 3), optional
            The periodic box vectors of the simulation in each frame, in nanometers.
            If not supplied, the vectors (1,0,0), (0,1,0) and (0,0,1) will
            be written for each frame.
        """
        if str(self.mode) != 'w':
            raise ValueError('write() is only available when the file is opened in mode="w"')

        # do typechecking, and then dispatch to the c level function
        xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False,
                          add_newaxis_on_deficient_ndim=True)
        n_frames = len(xyz)
        time = ensure_type(time, dtype=np.float32, ndim=1, name='time', can_be_none=True,
                           shape=(n_frames,), add_newaxis_on_deficient_ndim=True,
                           warn_on_cast=False)
        step = ensure_type(step, dtype=np.int32, ndim=1, name='step', can_be_none=True,
                           shape=(n_frames,), add_newaxis_on_deficient_ndim=True,
                            warn_on_cast=False)
        box = ensure_type(box, dtype=np.float32, ndim=3, name='box', can_be_none=True,
                          shape=(n_frames, 3, 3), add_newaxis_on_deficient_ndim=True,
                          warn_on_cast=False)
        if time is None:
            time = np.arange(0, n_frames, dtype=np.float32)
        if step is None:
            step = np.arange(0, n_frames, dtype=np.int32)
        if box is None:
            # make each box[i] be the identity matrix
            box = np.zeros((n_frames, 3, 3), dtype=np.float32)
            box[:,0,0] = np.ones(n_frames, dtype=np.float32)
            box[:,1,1] = np.ones(n_frames, dtype=np.float32)
            box[:,2,2] = np.ones(n_frames, dtype=np.float32)

        prec = 1000.0 * np.ones(n_frames, dtype=np.float32)

        if self.frame_counter == 0:
            self.n_atoms = xyz.shape[1]
        else:
            if not self.n_atoms == xyz.shape[1]:
                raise ValueError("This file has %d atoms, but you're now trying to write %d atoms" % (self.n_atoms, xyz.shape[1]))


        self._write(xyz, time, step, box, prec)


    def _write(self, np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz not None,
               np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time not None,
               np.ndarray[ndim=1, dtype=np.int32_t, mode='c'] step not None,
               np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box not None,
               np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] prec not None):

        cdef int n_frames = len(xyz)
        cdef int n_atoms = xyz.shape[1]
        cdef int status, i

        # all same shape
        assert n_frames == len(box) == len(step) == len(time) == len(prec)


        for i in range(n_frames):
            status = xdrlib.write_xtc(self.fh, n_atoms, step[i], time[i], <xdrlib.matrix>&box[i, 0, 0], <xdrlib.rvec*>&xyz[i, 0, 0], prec[i])
            if status != _EXDROK:
                raise RuntimeError('XTC write error: %s' % status)

        self.frame_counter += n_frames
        return status

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()
