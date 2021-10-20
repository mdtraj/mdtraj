# cython: c_string_type=str, c_string_encoding=ascii
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Teng Lin
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
import shutil
cimport numpy as np
np.import_array()
from mdtraj.utils import ensure_type, cast_indices, in_units_of
from mdtraj.utils.six import string_types
from mdtraj.formats.registry import FormatRegistry
from libc.stdlib cimport malloc, free
from dtrlib cimport molfile_timestep_t, molfile_timestep_metadata
from dtrlib cimport open_file_read, close_file_read, read_timestep2
from dtrlib cimport open_file_write, write_timestep, close_file_write
from dtrlib cimport read_timestep_metadata, read_times



##############################################################################
# Globals
##############################################################################

__all__ = ['DTRTrajectoryFile', 'load_dtr']

# codes that indicate status on return from library
cdef int _DTR_SUCCESS    = 0   # No problems
cdef int _DTR_EOF    = -1   # No problems

##############################################################################
# Code
##############################################################################


def _load_desmond_traj(filename, top=None, stride=None, atom_indices=None, frame=None):
    """
    """
    from mdtraj.core.trajectory import _parse_topology
    if top is None:
        raise ValueError('"top" argument is required for load_dtr')
    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_trr. '
            'you supplied %s' % type(filename))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)
    with DTRTrajectoryFile(str(filename)) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(topology, n_frames=n_frames, atom_indices=atom_indices, stride=stride)


@FormatRegistry.register_loader('.dtr')
def load_dtr(filename, top=None, stride=None, atom_indices=None, frame=None):
    """load_dtr(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load a dtr file from disk.

    The .dtr format is a cross-platform compressed binary trajectory format
    produced by DESMOND. Many different trajectory formats are generated
    from different versions of DESMOND. They usually stores:
        atomic coordinates
        velocities (or momentum that will be converted to velocity)
        box vectors
        time information
        other arbitrary data

    Below lists the supported trajectory format by various version of DESMOND:
        WRAPPED_V_2             : single precision, positions and velocities
        DBL_WRAPPED_V_2:        : double precision, positions and velocities
        WRAPPED_V_1             : single precision, positions and velocities
        DBL_WRAPPED_V_1         : double precision, positions and velocities
        POSN_MOMENTUM_V_1       : single precision, positions and momentum
        DBL_POSN_MOMENTUM_V_1   : double precision, positions and momentum
        ANTON_SFXP_V3           : trajectory generated on ANTON

    WRAPPED_V_2 is produced by DESMOND released by Schrodinger since 2010.

    Unlike other trajectory format, DESMOND trajectory splits frames into
    different files, and put them under a directory. Other auxiliary files
    containing meta data are also stored under the same directory. Among
    those files, there is one empty file named clickme.dtr, which is originally
    created to allow DESMOND trajectory to be loaded by VMD.

    Parameters
    ----------
    filename : path-like
        Path of dtr file that ends with clickme.dtr. Alternatively,
        directory name can also be accepted.
    top : {str, Trajectory, Topology}
        dtr format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
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
    >>> # loading clickme.dtr file under DESMOND trajectory directory
    >>> import mdtraj as md
    >>> traj = md.load_dtr('output_trj/clickme.dtr', top='topology.pdb')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>

    >>> # loading DESMOND trajectory directory directly
    >>> import mdtraj as md
    >>> traj = md.load_dtr('output_trj', top='topology.pdb')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.DTRTrajectoryFile :  Low level interface to dtr files
    """
    return _load_desmond_traj(filename, top=top, stride=stride, atom_indices=atom_indices, frame=frame)


@FormatRegistry.register_loader('.stk')
def load_stk(filename, top=None, stride=None, atom_indices=None, frame=None):

    """load_dtr(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load a stk file from disk.

    The .stk format is a text file where each line in the file points to a .dtr
    trajectory produced by DESMOND. It allow multiple trajectories to be
    concatenated. Note that, if chemical time of the listed trajectories
    are overlapped, the frames that are read in earlier with the same chemical
    time will be automatically discarded.

    Parameters
    ----------
    filename : path-like
        Path of stk file.
    top : {path-like, Trajectory, Topology}
        stk format does not contain topology information. Pass in either
        the path to a pdb file, a trajectory, or a topology to supply this
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
    >>> traj = md.load_stk('output.stk', top='topology.pdb')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.DTRTrajectoryFile :  Low level interface to dtr files
    """
    return _load_desmond_traj(filename, top=top, stride=stride, atom_indices=atom_indices, frame=frame)

cdef class DTRTrajectoryFile:
    """DTRTrajectoryFile(filename, mode='r', force_overwrite=True)

    Interface for reading and writing to a DESMOND dtr file.
    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    The conventional units in the dtr file are angstroms and degrees. The format
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
    >>> f = DTRTrajectoryFile('mytrajectory.dtr', 'r')
    >>> f.read(n_frames=1)  # read a single frame from the file
    >>> f.read()            # read all of the remaining frames
    >>> f.close()

    >>> # read all of the data with automatic closing of the file
    >>> with DTRTrajectoryFile('mytrajectory.dtr') as f:
    >>>    xyz, cell_lengths, cell_angles = f.read()

    >>> # write some xyz coordinates to a new file
    >>> with DTRTrajectoryFile('mytrajectory2.dtr. 'w') as f:
    >>>     f.write(np.random.randn(10,3,3))

    >>> # write frames one at a time
    >>> with DTRTrajectoryFile('mytrajectory2.dtr. 'w') as f:
    >>>     n_frames, n_atoms = 5, 10
    >>>     for i in range(n_frames):
    >>>         f.write(np.random.randn(n_atoms, 3))

    See Also
    --------
    mdtraj.load_dtr : High-level wrapper that returns a ``md.Trajectory``
    """

    # n_atoms and n_frames hold the number of atoms and the number of frames
    # in the file, as read off the header of the dtr file during read mode
    cdef int frame_counter, n_atoms, n_frames
    cdef void* fh
    cdef char* mode
    cdef char* filename
    cdef int is_open, _needs_write_initialization
    cdef molfile_timestep_t* timestep
    cdef molfile_timestep_metadata metadata

    cdef readonly char* distance_unit

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True):
        """Open a dtr Trajectory File
        """
        self.distance_unit = 'angstroms'
        self.is_open = False
        self.mode = mode

        if str(mode) == 'r':
            self.filename = filename
            self.fh = open_file_read(filename, "dtr", &self.n_atoms)
            if self.fh is NULL:
                raise IOError("Could not open file: %s" % filename)
            assert self.n_atoms > 0, 'dtr Corruption: n_atoms was not positive'
            read_timestep_metadata(self.fh, &self.metadata)
            self.n_frames = self.metadata.count
            self.frame_counter = 0
            self.is_open = True
        elif str(mode) == 'w':
            self.filename = filename
            self._needs_write_initialization = 1
            if os.path.exists(filename):
                if force_overwrite:
                    shutil.rmtree(filename)
                else:
                    raise IOError('"%s" already exists' % filename)
        else:
            raise ValueError("The only supported mode is 'r' ")

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
        """We don't actually want to open the dtr file during write mode
        until we know how many atoms the user wants to save. so this
        is delayed until the first call to write()
        """
        assert not self.is_open and self._needs_write_initialization

        self.n_atoms = n_atoms
        self.fh = open_file_write(self.filename, "dtr", self.n_atoms)
        if self.fh is NULL:
            raise IOError('There was an error opening the file: %s' % self.filename)
        self.is_open = True

        self._needs_write_initialization = False

    def get_times(self):
        """
        Reading chemical times associated with all frames.

        Returns
        -------
        times : np.ndarray, shape=(n_frames), dtype=float64
            The chemical time for each frame.
        """
        cdef np.ndarray[dtype=np.float64_t, ndim=1] times = np.zeros((self.n_frames), dtype=np.float64)
        count = read_times(self.fh, 0, self.n_frames, &times[0])
        if count != self.n_frames:
            raise IOError('Fail to read chemical time')
        return times


    def close(self):
        """Close the dtr file handle
        """
        if self.is_open and self.fh is not NULL:
            if str(self.mode) == 'r':
                close_file_read(self.fh)
            else:
                close_file_write(self.fh)
            self.is_open = False

        #self._needs_write_initialization = False


    def seek(self, offset, whence=0):
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
            absolute = self.n_frames + offset
        else:
            raise IOError('Invalid argument')

        if advance is not None:
            self.frame_counter += advance

        elif absolute is not None:
            self.frame_counter = absolute

        if self.frame_counter > self.n_frames:
            self.frame_counter = self.n_frames
        elif self.frame_counter < 0:
            self.frame_counter = 0

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
        This length is based on information in the header of the dtr
        file. It is possible for it to be off by 1 if, for instance,
        the client writing the file was killed between writing the last
        frame and updating the header information.
        """
        #pass
        if not self.is_open:
             raise ValueError('I/O operation on closed file')
        return self.n_frames

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """read_as_traj(topology, n_frames=None, stride=None, atom_indices=None)

        Read a trajectory from a DTR file

        Parameters
        ----------
        topology : Topology
            The system topology
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
        trajectory : Trajectory
            A trajectory object containing the loaded portion of the file.
        """
        from mdtraj.core.trajectory import Trajectory
        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        xyz, time, box_length, box_angle = self.read(stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)
        in_units_of(box_length, self.distance_unit, Trajectory._distance_unit, inplace=True)
        return Trajectory(xyz=xyz, topology=topology, time=time,
                          unitcell_lengths=box_length,
                          unitcell_angles=box_angle)

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """read(n_frames=None, stride=None, atom_indices=None)

        Read the data from a DTR file

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
            coordinates in the dtr file are stored in units of angstroms.
        times: np.array, shape=(n_frames), dtype=float64
            The chemical time in each frame.
        cell_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32
            The length of the box in each frame. `cell_lengths[i,0]` is the length
            of the A axis (in frame i), and `cell_lengths[i,1]` and
            `cell_lengths[i,2]` are the B and C axis respectively. By convention,
            the cell lengths in the dtr file are stored in units of angstroms.
        cell_angles : np.ndarray, shape=(n_frames, 3), dtype=float32
            Organized analogously to cell_lengths. Gives the alpha, beta and
            gamma angles respectively in entries `cell_angles[i,0]`,
            `cell_angles[i,1]`, `cell_angles[i,2]`. By convention, the cell
            angles in the dtr file are stored in units of degrees.
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

        if stride is None:
            _stride = 1
        else:
            _stride = stride

        _start = self.frame_counter
        _last = _start + _n_frames*_stride

        if _last > self.n_frames:
            _last = self.n_frames

        # # has trouble in  putting self.times into __init__
        # if not hasattr(self, 'times'):
        #     setattr(self, 'times', None)
        #
        # if not self.times:
        #     self.times = self.get_times()
        # TODO:
        times = self.get_times()
        times = times[_start:_last:_stride]

        _n_frames = len(times)

        if atom_indices is None:
            n_atoms_to_read = self.n_atoms
        elif isinstance(atom_indices, slice):
            n_atoms_to_read = len(np.arange(self.n_atoms)[atom_indices])
        else:
            if min(atom_indices) < 0:
                raise ValueError('atom_indices should be zero indexed. you gave an index less than zero')
            if max(atom_indices) >= self.n_atoms:
                raise ValueError('atom indices should be zero indexed. you gave an index bigger than the number of atoms')
            n_atoms_to_read = len(atom_indices)

        # allocate space to store the data that we are going to read off the disk
        # Desmond trajectory has different format, the storage could be very different
        # TODO: query the data type before the allocation
        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = np.zeros((_n_frames, n_atoms_to_read, 3), dtype=np.float32)
        cdef np.ndarray[dtype=np.float32_t, ndim=2] cell_lengths = np.zeros((_n_frames, 3), dtype=np.float32)
        cdef np.ndarray[dtype=np.float32_t, ndim=2] cell_angles = np.zeros((_n_frames, 3), dtype=np.float32)

        # only used if atom_indices is given
        cdef np.ndarray[dtype=np.float32_t, ndim=2] framebuffer = np.zeros((self.n_atoms, 3), dtype=np.float32)

        cdef int i, j
        cdef int status = _DTR_SUCCESS

        for j in range(_n_frames):
            i = j*_stride + _start

            if atom_indices is None:
                self.timestep.coords = &xyz[j,0,0]
                # set velocities to NULL, otherwise it will cause segmentation fault if the trajectory
                # happen to contain velocities
                self.timestep.velocities = NULL
                status = read_timestep2(self.fh, i, self.timestep)
            else:
                self.timestep.coords = &framebuffer[0,0]
                self.timestep.velocities = NULL
                status = read_timestep2(self.fh, i, self.timestep)
                xyz[j, :, :] = framebuffer[atom_indices, :]

            cell_lengths[j, 0] = self.timestep.A
            cell_lengths[j, 1] = self.timestep.B
            cell_lengths[j, 2] = self.timestep.C
            cell_angles[j, 0] = self.timestep.alpha
            cell_angles[j, 1] = self.timestep.beta
            cell_angles[j, 2] = self.timestep.gamma

            self.frame_counter += 1

            if status != _DTR_SUCCESS:
                raise IOError("Fail to read frame %d:"%i)

        if status == _DTR_SUCCESS:
            # if we're done either because of we read all of the n_frames
            # requested succcessfully, return
            return xyz, times, cell_lengths, cell_angles

        # if status == _DTR_EOF:
        #     # if we're doing because we reached a normal EOF (perhaps the)
        #     # user asked to read more frames than were in the file, we need
        #     # to truncate the return arrays -- we don't want to return them
        #     # a big stack of zeros.
        #     xyz = xyz[0:i]
        #     if cell_lengths is not None and cell_angles is not None:
        #         cell_lengths = cell_lengths[0:i]
        #         cell_angles = cell_angles[0:i]
        #
        #     return xyz, times, cell_lengths, cell_angles

        # If we got some other status, thats a "real" error.
        raise IOError("Error: %s", status)

    def write(self, xyz, cell_lengths=None, cell_angles=None, times=None):
        """write(xyz, cell_lengths=None, cell_angles=None, times=None)

        Write one or more frames of data to the dtr file

        Parameters
        ----------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms to write. By convention, the
            lengths should be in units of angstroms.
        cell_lengths : np.ndarray, shape=(n_frames, 3), dtype=float32
            The length of the periodic box in each frame, in each direction,
            `a`, `b`, `c`. By convention the lengths should be in units
            of angstroms.
        cell_angles : np.ndarray, shape=(n_frames, 3), dtype=float32
            Organized analogously to cell_lengths. Gives the alpha, beta and
            gamma angles respectively. By convention, the angles should be
            in units of degrees.
        times: np.ndarray, shape=(n_frames), dtype=float64
            The chemical time in each frame. It must be stored in ascending order.
        """


        if str(self.mode) != 'w':
            raise ValueError('write() is only available when the file is opened in mode="w"')
        if not self._needs_write_initialization and not self.is_open:
            raise IOError("file is not open")


        # do typechecking, and then dispatch to the c level function
        xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False,
                          add_newaxis_on_deficient_ndim=True)
        n_frames = len(xyz)

        # they must be both present or both absent
        if cell_lengths is None or cell_angles is None or times is None:
            raise ValueError('cell_lengths, cell_angles and times must be given')

        # check if times is in ascending order
        sorted = times.argsort()
        if not np.all(sorted[1:]>sorted[:-1]):
            raise ValueError("times must be in ascending order")


        cell_lengths = ensure_type(cell_lengths, dtype=np.float32, ndim=2, name='cell_lengths',
                                   can_be_none=False, shape=(n_frames, 3), add_newaxis_on_deficient_ndim=True,
                                   warn_on_cast=False)
        cell_angles = ensure_type(cell_angles, dtype=np.float32, ndim=2, name='cell_angles',
                                  can_be_none=False, shape=(n_frames, 3), add_newaxis_on_deficient_ndim=True,
                                  warn_on_cast=False)

        times = ensure_type(times, dtype=np.float64, ndim=1, name='times', can_be_none=False,
                          add_newaxis_on_deficient_ndim=True)

        if self._needs_write_initialization:
            self._initialize_write(xyz.shape[1])
        else:
            if not self.n_atoms == xyz.shape[1]:
                raise ValueError('Number of atoms doesnt match')
            if cell_lengths is None and cell_angles:
                raise ValueError("The file that you're saving to expects each frame "
                    "to contain unitcell information, but you did not supply it.")
            if times is None:
                raise ValueError("The file that you're saving to was created without "
                    "time information.")

        self._write(xyz, cell_lengths, cell_angles, times)


    cdef _write(self, np.ndarray[np.float32_t, ndim=3, mode="c"] xyz,
                np.ndarray[np.float32_t, ndim=2, mode="c"] cell_lengths,
                np.ndarray[np.float32_t, ndim=2, mode="c"] cell_angles,
                np.ndarray[np.float64_t, ndim=1, mode="c"] times):
        if str(self.mode) != 'w':
            raise ValueError('_write() is only available when the file is opened in mode="w"')


        cdef int i, status
        cdef int n_frames = len(xyz)

        for i in range(n_frames):
            self.timestep.coords = &xyz[i, 0, 0]
            self.timestep.velocities = NULL
            if cell_angles is not None and cell_lengths is not None:
                self.timestep.A = cell_lengths[i, 0]
                self.timestep.B = cell_lengths[i, 1]
                self.timestep.C = cell_lengths[i, 2]
                self.timestep.alpha = cell_angles[i, 0]
                self.timestep.beta  = cell_angles[i, 1]
                self.timestep.gamma = cell_angles[i, 2]
                self.timestep.physical_time = times[i]

            # when the dcd handle is opened during initialize_write,
            # it passes the flag for whether the unitcell information is written
            # to disk or not. So when we don't have unitcell information being
            # written, the unitcell fields (A,B,C,alpha,beta,gamma) are ignored
            # during write_timestep()

            status = write_timestep(self.fh, self.timestep)

            if status != _DTR_SUCCESS:
                raise IOError("DTR Error: %s" %status)



FormatRegistry.register_fileobject('.dtr')(DTRTrajectoryFile)
