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


##############################################################################
# Imports
##############################################################################

import cython
cimport cython
import os
import numpy as np
cimport numpy as np
np.import_array()
from mdtraj.utils.arrays import ensure_type
from libc.stdlib cimport malloc, free
from dcdlib cimport molfile_timestep_t, dcdhandle
from dcdlib cimport open_dcd_read, close_file_read, read_next_timestep
from dcdlib cimport open_dcd_write, close_file_write, write_timestep
from dcdlib cimport dcd_nsets


##############################################################################
# Globals
##############################################################################

# codes that indicate status on return from library
cdef int _DCD_SUCCESS    = 0   # No problems
cdef int _DCD_EOF    = -1   # No problems


cdef ERROR_MESSAGES = {
    -1: 'Normal EOF',
    -2: 'DCD file does not exist',
    -3: 'Open of DCD file failed',
    -4: 'Read call on DCD file failed',
    -5: 'Premature EOF found in DCD file',
    -6: 'Format of DCD file is wrong',
    -7: 'Output file already exists',
    -8: 'Malloc failed',
    -9: 'Write call on DCD file failed',
}

##############################################################################
# Classes
##############################################################################


cdef class DCDTrajectoryFile:
    """DCDTrajectoryFile(filename, mode='r', force_overwrite=True)

    Interface for reading and writing to a CHARMM/NAMD DCD file.
    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    The conventional units in the DCD file are angstroms and degrees. The format
    supports saving coordinates and unit cell parameters (lengths and angles)

    Parameters
    ----------
    filename : string
        Path to the file to open
    mode : {'r', 'w'}
        Mode in which to open the file. 'r' is for reading, and 'w' is for writing.
    force_overwrite : bool
        In mode='w', how do you want to behave if a file by the name of `filename`
        already exists? if `force_overwrite=True`, it will be overwritten.

    Examples
    --------
    >>> # read a single frame, and then read the remaining frames
    >>> f = DCDTrajectoryFile('mytrajectory.dcd', 'r')
    >>> f.read(n_frames=1)  # read a single frame from the file
    >>> xyzf.read()            # read all of the remaining frames
    >>> f.close()

    >>> # read all of the data with automatic closing of the file
    >>> with DCDTrajectoryFile('mytrajectory.dcd') as f:
    >>>    xyz, cell_lengths, cell_angles = f.read()

    >>> # write some xyz coordinates to a new file
    >>> with DCDTrajectoryFile('mytrajectory2.dcd. 'w') as f:
    >>>     f.write(np.random.randn(10,3,3))

    >>> # write frames one at a time
    >>> with DCDTrajectoryFile('mytrajectory2.dcd. 'w') as f:
    >>>     n_frames, n_atoms = 5, 10
    >>>     for i in range(n_frames):
    >>>         f.write(np.random.randn(n_atoms, 3))

    See Also
    --------
    mdtraj.load_dcd : High-level wrapper that returns a ``md.Trajectory``
    """

    # n_atoms and n_frames hold the number of atoms and the number of frames
    # in the file, as read off the header of the DCD file during read mode
    cdef int frame_counter, n_atoms, n_frames
    cdef dcdhandle* fh
    cdef char* mode
    cdef char* filename
    cdef int is_open, _needs_write_initialization
    cdef molfile_timestep_t* timestep
    cdef readonly char* distance_unit

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True):
        """Open a DCD Trajectory File
        """
        self.distance_unit = 'angstroms'
        self.is_open = False
        self.mode = mode

        if str(mode) == 'r':
            self.filename = filename
            self.fh = open_dcd_read(filename, "dcd", &self.n_atoms, &self.n_frames)
            assert self.n_atoms > 0, 'DCD Corruption: n_atoms was not positive'
            assert self.n_frames >= 0, 'DCD corruption: n_frames < 0'
            # we're at the beginning of the file now
            self.frame_counter = 0
            self.is_open = True
            if self.fh is NULL:
                raise IOError('There was an error opening the file: %s' % filename)
        elif str(mode) == 'w':
            self.filename = filename
            self._needs_write_initialization = 1
            if not force_overwrite and os.path.exists(filename):
                raise IOError('"%s" already exists' % filename)
        else:
            raise ValueError("most must be one of ['r', 'w']")

        # alloc the molfile_timestep, which is the struct that the library is
        # going to deposit its data into each timestep
        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.timestep is NULL:
            raise MemoryError('There was an error allocating memory')

    def __dealloc__(self):
        # free whatever we malloced
        free(self.timestep)
        self.close()

    def _initialize_write(self, int n_atoms):
        """We don't actually want to open the dcd file during write mode
        until we know how many atoms the user wants to save. so this
        is delayed until the first call to write()
        """
        assert not self.is_open and self._needs_write_initialization

        self.n_atoms = n_atoms
        self.fh = open_dcd_write(self.filename, "dcd", self.n_atoms)
        if self.fh is NULL:
            raise IOError('There was an error opening the file: %s' % self.filename)
        self.is_open = True

        self._needs_write_initialization = False

    def close(self):
        "Close the DCD file handle"
        if self.is_open and self.fh is not NULL:
            if str(self.mode) == 'r':
                close_file_read(self.fh)
            else:
                close_file_write(self.fh)
            self.is_open = False

        self._needs_write_initialization = False

    def seek(self, offset, whence=0):
        """Move to a new file position

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
        cdef int i, status
        if str(self.mode) != 'r':
            raise NotImplementedError("seek is only supported in mode='r'")

        advance, absolute = None, None
        if whence == 0 and offset >= 0:
            if offset >= self.frame_counter:
                advance = offset - self.frame_counter
            else:
                absolute = offset
        elif whence == 1 and offset >= 0:
            advance = offset
        elif whence == 1 and offset < 0:
            absolute = offset + self.frame_counter
        elif whence == 2 and offset <= 0:
            raise NotImplementedError('offsets from the end are not supported yet')
        else:
            raise IOError('Invalid argument')

        if advance is not None:
            self.frame_counter += advance
            for i in range(advance):
                status = read_next_timestep(self.fh, self.n_atoms, NULL)
        elif absolute is not None:
            close_file_read(self.fh)
            self.fh = open_dcd_read(self.filename, "dcd", &self.n_atoms, &self.n_frames)
            for i in range(absolute):
                status = read_next_timestep(self.fh, self.n_atoms, NULL)
            self.frame_counter = absolute

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        return int(self.frame_counter)

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

    def __len__(self):
        """Number of frames in the file

        Notes
        -----
        This length is based on information in the header of the DCD
        file. It is possible for it to be off by 1 if, for instance,
        the client writing the file was killed between writing the last
        frame and updating the header information.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        return dcd_nsets(self.fh)

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """read(n_frames=None, stride=None, atom_indices=None)

        Read the data from a DCD file

        Parameters
        ----------
        n_frames : int, optional
            If positive, then read only the next `n_frames` frames. Otherwise read all
            of the frames in the file.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it required
            an extra copy, but will save memory.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=float32
            The xyz coordinates of each atom in each frame. By convention, the
            coordinates in the dcd file are stored in units of angstroms.
        cell_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32
            The length of the box in each frame. `cell_lengths[i,0]` is the length
            of the A axis (in frame i), and `cell_lengths[i,1]` and
            `cell_lengths[i,2]` are the B and C axis respectively. By convention,
            the cell lengths in the dcd file are stored in units of angstroms.
        cell_angles : np.ndarray, shape=(n_frames, 3), dtype=float32
            Organized analogously to cell_lengths. Gives the alpha, beta and
            gamma angles respectively in entries `cell_angles[i,0]`,
            `cell_angles[i,1]`, `cell_angles[i,2]`. By convention, the cell
            angles in the dcd file are stored in units of degrees.
        """
        if str(self.mode) != 'r':
            raise ValueError('read() is only available when the file is opened in mode="r"')
        if not self.is_open:
            raise IOError("file is not open")

        cdef int _n_frames, n_atoms_to_read, _stride
        if n_frames is None:
            # if the user specifies n_frames=None, they want to read to the
            # end of the file
            _n_frames = self.n_frames - self.frame_counter
        else:
            _n_frames = int(n_frames)

        if atom_indices is None:
            n_atoms_to_read = self.n_atoms
        elif isinstance(atom_indices, slice):
            n_atoms_to_read = len(np.arange(self.n_atoms)[atom_indices])
        else:
            if min(atom_indices) < 0:
                raise ValueError('atom_indices should be zero indexed. you gave an index less than zerp')
            if max(atom_indices) >= self.n_atoms:
                raise ValueError('atom indices should be zero indexed. you gave an index bigger than the number of atoms')
            n_atoms_to_read = len(atom_indices)

        if stride is None:
            _stride = 1
        else:
            _stride = stride

        # malloc space to put the data that we're going to read off the disk
        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = np.zeros((_n_frames, n_atoms_to_read, 3), dtype=np.float32)
        cdef np.ndarray[dtype=np.float32_t, ndim=2] cell_lengths = np.zeros((_n_frames, 3), dtype=np.float32)
        cdef np.ndarray[dtype=np.float32_t, ndim=2] cell_angles = np.zeros((_n_frames, 3), dtype=np.float32)

        # only used if atom_indices is given
        cdef np.ndarray[dtype=np.float32_t, ndim=2] framebuffer = np.zeros((self.n_atoms, 3), dtype=np.float32)

        cdef int i, j
        cdef int status = _DCD_SUCCESS

        for i in range(_n_frames):
            if atom_indices is None:
                self.timestep.coords = &xyz[i,0,0]
                status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
            else:
                self.timestep.coords = &framebuffer[0,0]
                status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
                xyz[i, :, :] = framebuffer[atom_indices, :]

            self.frame_counter += 1
            cell_lengths[i, 0] = self.timestep.A
            cell_lengths[i, 1] = self.timestep.B
            cell_lengths[i, 2] = self.timestep.C
            cell_angles[i, 0] = self.timestep.alpha
            cell_angles[i, 1] = self.timestep.beta
            cell_angles[i, 2] = self.timestep.gamma

            if status != _DCD_SUCCESS:
                # if the frame was not successfully read, then we're done
                break

            for j in range(_stride - 1):
                status = read_next_timestep(self.fh, self.n_atoms, NULL)

        if status == _DCD_SUCCESS:
            # if we're done either because of we read all of the n_frames
            # requested succcessfully, return
            return xyz, cell_lengths, cell_angles

        if status == _DCD_EOF:
            # if we're doing because we reached a normal EOF (perhaps the)
            # user asked to read more frames than were in the file, we need
            # to truncate the return arrays -- we don't want to return them
            # a big stack of zeros.
            xyz = xyz[0:i]
            cell_lengths = cell_lengths[0:i]
            cell_angles = cell_angles[0:i]

            return xyz, cell_lengths, cell_angles

        # If we got some other status, thats a "real" error.
        raise IOError("Error: %s", ERROR_MESSAGES(status))

    def write(self, xyz, cell_lengths=None, cell_angles=None):
        """write(xyz, cell_lengths=None, cell_angles=None)

        Write one or more frames of data to the DCD file

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. By convention, the
            lengths should be in units of angstroms.
        cell_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
            The length of the periodic box in each frame, in each direction,
            `a`, `b`, `c`. By convention the lengths should be in units
            of angstroms.
        cell_angles : np.ndarray, shape=(n_frames, 3), dtype=float32, optional
            Organized analogously to cell_lengths. Gives the alpha, beta and
            gamma angles respectively. By convention, the angles should be
            in units of degrees.
        """
        if str(self.mode) != 'w':
            raise ValueError('write() is only available when the file is opened in mode="w"')
        if not self._needs_write_initialization and not self.is_open:
            raise IOError("file is not open")

        # do typechecking, and then dispatch to the c level function
        xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False,
                          add_newaxis_on_deficient_ndim=True)
        n_frames = len(xyz)

        cell_lengths = ensure_type(cell_lengths, dtype=np.float32, ndim=2, name='cell_lengths',
                                   can_be_none=True, shape=(n_frames, 3), add_newaxis_on_deficient_ndim=True,
                                   warn_on_cast=False)
        if cell_lengths is None:
            cell_lengths = np.ones((n_frames, 3), dtype=np.float32)

        cell_angles = ensure_type(cell_angles, dtype=np.float32, ndim=2, name='cell_angles',
                                  can_be_none=True, shape=(n_frames, 3), add_newaxis_on_deficient_ndim=True,
                                  warn_on_cast=False)
        if cell_angles is None:
            cell_angles = 90.0 * np.ones((n_frames, 3), dtype=np.float32)


        if self._needs_write_initialization:
            self._initialize_write(xyz.shape[1])
        else:
            if not self.n_atoms == xyz.shape[1]:
                raise ValueError('Number of atoms doesnt match')
        self._write(xyz, cell_lengths, cell_angles)

    cdef _write(self, np.ndarray[np.float32_t, ndim=3, mode="c"] xyz,
                np.ndarray[np.float32_t, ndim=2, mode="c"] cell_lengths,
                np.ndarray[np.float32_t, ndim=2, mode="c"] cell_angles):
        if str(self.mode) != 'w':
            raise ValueError('_write() is only available when the file is opened in mode="w"')


        cdef int i, status
        cdef int n_frames = len(xyz)

        for i in range(n_frames):
            self.timestep.coords = &xyz[i, 0, 0]
            self.timestep.A = cell_lengths[i, 0]
            self.timestep.B = cell_lengths[i, 1]
            self.timestep.C = cell_lengths[i, 2]
            self.timestep.alpha = cell_angles[i, 0]
            self.timestep.beta  = cell_angles[i, 1]
            self.timestep.gamma = cell_angles[i, 2]

            status = write_timestep(self.fh, self.timestep)

            if status != _DCD_SUCCESS:
                raise IOError("DCD Error: %s" % ERROR_MESSAGES(status))
