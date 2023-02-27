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
from __future__ import print_function

import os
import warnings
import xdrlib
import numpy as np
cimport numpy as np
np.import_array()

from mdtraj.utils import ensure_type, cast_indices, in_units_of
from mdtraj.utils.six import string_types
from mdtraj.formats.registry import FormatRegistry

cimport xdrlib

from libc.stdio cimport SEEK_SET, SEEK_CUR
from libc.math cimport ceil

ctypedef np.npy_int64   int64_t


__all__ = ['load_xtc', 'XTCTrajectoryFile']


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

# Note: following constant depend on int is a 32bit integer!

# numpy variable types include the specific numpy of bytes of each, but the c
# variables in our interface file don't. this could get bad if we're on a wierd
# machine, so lets make sure first
if sizeof(int) != sizeof(np.int32_t):
    raise RuntimeError('Integers on your compiler are not 32 bits. This is not good.')
if sizeof(float) != sizeof(np.float32_t):
    raise RuntimeError('Floats on your compiler are not 32 bits. This is not good')

# constants for short (<= 10 atoms) XTC files:
cdef int DIM = 3
# header fields (before coords):
# 1. magic (int)
# 2. natoms (int)
# 3. step (int)
# 4. time (float)
cdef int XTC_HDR_SIZE = 3*sizeof(np.int32_t) + sizeof(np.float32_t)
cdef int XTC_SHORT_HEADER_SIZE = XTC_HDR_SIZE + DIM**2 * sizeof(np.float32_t) + 4
cdef int XTC_SHORT_BYTES_PER_ATOM = DIM*sizeof(np.float32_t)

# constant for 'regular' XTCs (> 10 atoms):
# 1. magic(int)
# 2. natoms (int)
# 3. step (int)
# 4. time (float)
# 5. DIM*DIM_box_vecs (DIM*DIM floats)
# 6. natoms (int)
# 7. prec (float)
# 8. DIM_min_xyz (int[3])
# 9. DIM_max_xyz (int[3])
# 10. smallidx (int)
cdef int XTC_HEADER_SIZE = 11*sizeof(np.int32_t) + 2*sizeof(np.float32_t) + DIM**2 * sizeof(np.float32_t)


###############################################################################
# Code
###############################################################################

@FormatRegistry.register_loader('.xtc')
def load_xtc(filename, top=None, stride=None, atom_indices=None, frame=None):
    """load_xtc(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load a Gromacs XTC file from disk.

    The .xtc format is a cross-platform compressed binary trajectory format
    produced by the gromacs software that stores atomic coordinates, box
    vectors, and time information. It is lossy (storing coordinates to about
    1e-3 A) and extremely space-efficient.

    Parameters
    ----------
    filename : path-like
        Filename (string) of xtc trajectory.
    top : {str, Trajectory, Topology}
        The XTC format does not contain topology information. Pass in either the
        path to a RCSB PDB file, a trajectory, or a topology to supply this
        information.
    stride : int, default=None
        Only read every stride-th frame
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. This may be slightly slower than the standard read because it
        requires an extra copy, but will save memory.
    frame : int, optional
        Use this option to load only a single frame from a trajectory on disk.
        If frame is None, the default, the entire trajectory will be loaded.
        If supplied, ``stride`` will be ignored.

    Examples
    --------
    >>> import mdtraj as md
    >>> traj = md.load_xtc('output.xtc', top='topology.pdb')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.XTCTrajectoryFile :  Low level interface to XTC files
    """
    from mdtraj.core.trajectory import _parse_topology
    if top is None:
        raise ValueError('"top" argument is required for load_xtc')

    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_xtc. '
                        'you supplied %s' % type(filename))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with XTCTrajectoryFile(str(filename), 'r') as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(topology, n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)


cdef class XTCTrajectoryFile(object):
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
    cdef str filename
    cdef int n_atoms          # number of atoms in the file
    cdef int64_t n_frames # number of frames in the file, cached
    cdef int64_t frame_counter    # current position in the file, in read mode
    cdef char is_open          # is the file handle currently open?
    cdef int64_t approx_n_frames  # appriximate number of frames in the file, as guessed based on its size
    cdef char* mode           # mode in which the file is open, either 'r' or 'w'
    cdef int min_chunk_size
    cdef float chunk_size_multiplier
    cdef char with_unitcell    # used in mode='w' to know if we're writing unitcells or nor
    cdef readonly char* distance_unit
    cdef np.ndarray _offsets

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True, **kwargs):
        """Open a GROMACS XTC file for reading/writing.
        """
        self.distance_unit = 'nanometers'
        self.is_open = False
        self.frame_counter = 0
        self.n_frames = -1  # means unknown
        self.filename = filename
        self._offsets = None

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
            if not force_overwrite and os.path.exists(filename):
                raise IOError('"%s" already exists' % filename)
            self.fh = xdrlib.xdrfile_open(str(filename), 'w')
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

        return max(approx_n_frames, 1)

    def __dealloc__(self):
        self.close()

    def close(self):
        "Close the XTC file handle"
        if self.is_open:
            xdrlib.xdrfile_close(self.fh)
            self.is_open = False

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """read_as_traj(topology, n_frames=None, stride=None, atom_indices=None)

        Read a trajectory from an XTC file

        Parameters
        ----------
        topology : Topology
            The system topology
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
        trajectory : Trajectory
            A trajectory object containing the loaded portion of the file.

        See Also
        --------
        read : Returns the raw data from the file
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        xyz, time, step, box = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(box, self.distance_unit, Trajectory._distance_unit, inplace=True)

        trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
        trajectory.unitcell_vectors = box
        return trajectory

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
        status: int
            either _EXDROK or _EXDRENDOFFILE

        See Also
        --------
        read_as_traj : Returns a Trajectory object
        """
        if not str(self.mode) == 'r':
            raise ValueError('read() is only available when file is opened in mode="r"')
        if not self.is_open:
            raise IOError('file must be open to read from it.')
        stride = int(stride) if stride is not None else 1
        if n_frames is not None:
            # if they supply the number of frames they want, that's easy
            if not int(n_frames) == n_frames:
                raise ValueError('n_frames must be an int, you supplied "%s"' % n_frames)
            xyz, time, step, box, status = self._read(int(n_frames), atom_indices, stride)
            if np.all(np.logical_and(box < 1e-10, box > -1e-10)):
                box = None
            return xyz, time, step, box

        # if they want ALL of the remaining frames, we need to guess at the
        # chunk size, and then check the exit status to make sure we're really
        # at the EOF
        all_xyz, all_time, all_step, all_box = [], [], [], []
        status = _EXDROK
        while status == _EXDROK:
            # guess the size of the chunk to read, based on how many frames we
            # think are in the file and how many we've currently read
            chunk = max(abs(int((self.approx_n_frames - self.frame_counter) * self.chunk_size_multiplier / stride)),
                        self.min_chunk_size)
            xyz, time, step, box, status = self._read(chunk, atom_indices, stride)

            all_xyz.append(xyz)
            all_time.append(time)
            all_step.append(step)
            all_box.append(box)

        if len(all_xyz) == 0:
            return np.empty(0), np.empty(0), np.empty(0), np.empty(0)
        all_xyz = np.concatenate(all_xyz)
        all_time = np.concatenate(all_time)
        all_step = np.concatenate(all_step)
        all_box =  np.concatenate(all_box)
        if np.all(np.logical_and(all_box < 1e-10, all_box > -1e-10)):
            all_box = None
        return all_xyz, all_time, all_step, all_box

    def _read(self, int64_t n_frames, atom_indices, stride):
        """Read a specified number of XTC frames from the buffer"""
        cdef int64_t n_read_frames = 0
        cdef int status = _EXDROK
        cdef int status_seek = _EXDROK
        cdef unsigned int n_atoms_to_read
        cdef char efficient_striding = stride > 1 and self._offsets is not None

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
            np.empty(n_frames, dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.int32_t, mode='c'] step = \
            np.empty(n_frames, dtype=np.int32)
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box = \
            np.empty((n_frames, 3, 3), dtype=np.float32)
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] prec = \
            np.empty(n_frames, dtype=np.float32)

        # only used if atom_indices is given
        cdef np.ndarray[dtype=np.float32_t, ndim=2] framebuffer
        if atom_indices is not None:
            framebuffer = np.empty((self.n_atoms, 3), dtype=np.float32)

        # striding dummy, only used if efficient_striding is false or at the end of the file.
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz_stride
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] time_stride
        cdef np.ndarray[ndim=1, dtype=np.int32_t, mode='c'] step_stride
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box_stride
        cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] prec_stride

        if stride > 1 and not efficient_striding:
            xyz_stride = np.empty((1, self.n_atoms, 3), dtype=np.float32)
            time_stride = np.empty(1, dtype=np.float32)
            step_stride = np.empty(1, dtype=np.int32)
            box_stride = np.empty((1, 3, 3), dtype=np.float32)
            prec_stride = np.empty(1, dtype=np.float32)

        while (n_read_frames < n_frames) and (status != _EXDRENDOFFILE):
            if atom_indices is None:
                status = xdrlib.read_xtc(self.fh, self.n_atoms, <int*> &step[n_read_frames],
                                         &time[n_read_frames], <xdrlib.matrix>&box[n_read_frames,0,0], <xdrlib.rvec*>&xyz[n_read_frames,0,0], &prec[n_read_frames])
            else:
                status = xdrlib.read_xtc(self.fh, self.n_atoms, <int*> &step[n_read_frames],
                                         &time[n_read_frames], <xdrlib.matrix>&box[n_read_frames,0,0], <xdrlib.rvec*>&framebuffer[0,0], &prec[n_read_frames])
                xyz[n_read_frames, :, :] = framebuffer[atom_indices, :]

            if status != _EXDRENDOFFILE and status != _EXDROK:
                raise RuntimeError('XTC read error: %s' % _EXDR_ERROR_MESSAGES.get(status, 'unknown'))
            else:
                # we have successfully read a frame!
                n_read_frames += 1

                if stride > 1:
                    # Can we seek within bounds?
                    if efficient_striding:
                        if self.frame_counter + stride < len(self):
                            self.seek(stride, whence=1)
                        else:
                            # deliberately dont set eof status as last frame is valid!
                            status = _EXDRENDOFFILE
                            n_read_frames += 1
                    else:
                        for _ in range(stride - 1):
                            seek_status = xdrlib.read_xtc(self.fh, self.n_atoms, <int*> &step_stride[0],
                                                          &time_stride[0], <xdrlib.matrix> &box_stride[0,0,0],
                                                          <xdrlib.rvec*>&xyz_stride[0,0,0], &prec_stride[0])
                            if seek_status != _EXDROK:
                                break

        if status == _EXDRENDOFFILE:
            xyz = xyz[:n_read_frames-1]
            box = box[:n_read_frames-1]
            time = time[:n_read_frames-1]
            step = step[:n_read_frames-1]

        # if we are using seek, the frame_counter already points to the right absolute position,
        # otherwise we increment the counter relatively
        if not efficient_striding:
            self.frame_counter += len(xyz)

        return xyz, time, step, box, status

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
                          add_newaxis_on_deficient_ndim=True, warn_on_cast=False)
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

        if self.frame_counter == 0:
            self.n_atoms = xyz.shape[1]
            self.with_unitcell = (box is not None)
        else:
            if not self.n_atoms == xyz.shape[1]:
                raise ValueError("This file has %d atoms, but you're now "
                    "trying to write %d atoms" % (self.n_atoms, xyz.shape[1]))
            if self.with_unitcell and (box is None):
                raise ValueError("The file that you're saving to expects each frame "
                    "to contain unitcell information, but you did not supply it.")
            if (not self.with_unitcell) and (box is not None):
                raise ValueError("The file that you're saving to was created without "
                    "unitcell information.")

        if time is None:
            time = np.arange(0, n_frames, dtype=np.float32)
        if step is None:
            step = np.arange(0, n_frames, dtype=np.int32)
        if box is None:
            # make each box[i] be the all zeros, which indicates the lack of
            # a unitcell
            box = np.zeros((n_frames, 3, 3), dtype=np.float32)

        prec = 1000.0 * np.ones(n_frames, dtype=np.float32)
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

    def seek(self, int64_t offset, int whence=0):
        """seek(offset, whence=0)

        Move to a new file position

        Parameters
        ----------
        offset : int
            A number of frames.
        whence : {0, 1, 2}
            0: offset from start of file, offset should be >=0.
            1: move relative to the current position, positive or negative
            2: move relative to the end of file, offset should be <= 0.
            Seeking beyond the end of a file is not supported
        """
        cdef int status
        cdef int64_t pos, absolute

        if str(self.mode) != 'r':
            raise NotImplementedError('seek() only available in mode="r" currently')
        if whence == 0 and offset >= 0:
            absolute = offset
        elif whence == 1:
            absolute = offset + self.frame_counter
        elif whence == 2 and offset <= 0:
            raise NotImplementedError('offsets from the end are not supported yet')
        else:
            raise IOError('Invalid argument')

        if absolute < 0 or absolute >= len(self.offsets):
            raise IOError('XTC Seek out of bounds: given absolute position: {}'.format(absolute))

        pos = self.offsets[absolute]
        status = xdrlib.xdr_seek(self.fh, pos, SEEK_SET)
        if status != _EXDROK:
            raise RuntimeError('XTC seek error: %s' % status)
        self.frame_counter = absolute

    def _calc_len_and_offsets(self):
        cdef int byte_offset, status
        cdef int64_t n_frames, filesize
        cdef np.ndarray[dtype=int64_t, ndim=1] offsets

        # restore old pos when done or in case of error
        cdef int64_t old_pos = xdrlib.xdr_tell(self.fh)

        if self.n_atoms <= 9:
            filesize = os.stat(self.filename).st_size
            byte_offset = XTC_SHORT_HEADER_SIZE + XTC_SHORT_BYTES_PER_ATOM*self.n_atoms
            assert filesize % byte_offset == 0, ("filesize(%i) not divideable"
                                                 " by bytes per frames(%i)"
                                                 % (filesize, byte_offset))
            n_frames = filesize // byte_offset
            offsets = np.fromiter((i*byte_offset for i in range(n_frames)),
                                  dtype=np.int64, count=n_frames)
        else:
            offsets = np.empty(self.approx_n_frames, dtype=np.int64)
            assert len(offsets) >= 1

            try:
                # skip header
                if xdrlib.xdr_seek(self.fh, XTC_HEADER_SIZE, SEEK_SET) != 0:
                    raise RuntimeError('could not skip header for file ' + self.filename)

                # init first byte_offset
                status = xdrlib.xdrfile_read_int(&byte_offset, 1, self.fh)
                if status == 0:
                    raise RuntimeError("error reading from first frame: %i" % status)
                byte_offset += 3 - ((byte_offset + 3) % 0x04)

                n_frames = 1
                offsets[0] = 0
                resize = np.resize
                while True:
                    # relative seek
                    status = xdrlib.xdr_seek(self.fh, byte_offset + XTC_HEADER_SIZE, SEEK_CUR)
                    if status != 0:
                        offset = byte_offset + XTC_HEADER_SIZE
                        last_pos = xdrlib.xdr_tell(self.fh)
                        raise RuntimeError("error during seek: status code "
                                           "fseek=%s; offset=%s, last_pos=%s"
                                           % (status, offset, last_pos))

                    # return value == # ints read, so we're finished
                    if xdrlib.xdrfile_read_int(&byte_offset, 1, self.fh) == 0:
                        break

                    if n_frames == len(offsets):
                        new_len = int(ceil(len(offsets)*1.2))
                        offsets = resize(offsets, new_len)

                    offsets[n_frames] = xdrlib.xdr_tell(self.fh) - 4 - XTC_HEADER_SIZE
                    n_frames += 1

                    # truncate byte offset to next 32-bit boundary
                    byte_offset += 3 - ((byte_offset + 3) % 0x04)
            finally:
                xdrlib.xdr_seek(self.fh, old_pos, SEEK_SET)

        return n_frames, offsets[:n_frames]

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        if str(self.mode) != 'r':
            raise NotImplementedError('tell() only available in mode="r" currently')
        return int(self.frame_counter)

    @property
    def offsets(self):
        "get byte offsets from current xtc file"
        if self._offsets is None:
            self.n_frames, self._offsets = self._calc_len_and_offsets()
        return self._offsets

    @offsets.setter
    def offsets(self, offsets):
        "set frame offsets"
        self._offsets = offsets

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

    def __len__(self):
        "Number of frames in the file"
        if str(self.mode) != 'r':
            raise NotImplementedError('len() only available in mode="r" currently')
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.n_frames == -1:
            self.offsets # invokes _calc_len_and_offsets
        return int(self.n_frames)

    def flush(self):
        if str(self.mode) != 'w':
            raise RuntimeError('could not flush file opened only for reading')
        cdef int status = xdrlib.xdr_flush(self.fh)
        if status != 0:
            raise IOError('could not flush xtc file. Status: %s' % status)


FormatRegistry.register_fileobject('.xtc')(XTCTrajectoryFile)
