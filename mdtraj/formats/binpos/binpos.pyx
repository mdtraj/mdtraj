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

import cython
import warnings
cimport cython
from libc.stdio cimport SEEK_SET, SEEK_CUR, SEEK_END
import os
import numpy as np
cimport numpy as np
np.import_array()
from mdtraj.utils import ensure_type, cast_indices, in_units_of
from mdtraj.utils.six import string_types
from mdtraj.formats.registry import FormatRegistry
from libc.stdlib cimport malloc, free
from binposlib cimport molfile_timestep_t
from binposlib cimport seek_timestep, tell_timestep;
from binposlib cimport open_binpos_read, close_file_read, read_next_timestep
from binposlib cimport open_binpos_write, close_file_write, write_timestep

###############################################################################
# Globals
###############################################################################

__all__ = ['BINPOSTrajectoryFile', 'load_binpos']
# codes that indicate status on return from library
cdef int _BINPOS_SUCESS = 0  # regular exit code
cdef int _BINPOS_EOF = -1  # end of file (or error)

###############################################################################
# Classes
###############################################################################

@FormatRegistry.register_loader('.binpos')
def load_binpos(filename, top=None, stride=None, atom_indices=None, frame=None):
    """load_binpos(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load an AMBER .binpos file from disk.

    The .binpos format is a cross-platform binary trajectory format produced by
    AMBER software. It stores only the atomic coordinates, and does *not* store
    any unit cell informations. Its use is discouraged.

    Parameters
    ----------
    filename : path-like
        Path of AMBER binpos file.
    top : {path-like, Trajectory, Topology}
        The BINPOS format does not contain topology information. Pass in either
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
    >>> import mdtraj as md                                        # doctest: +SKIP
    >>> traj = md.load_binpos('output.binpos', top='topology.pdb') # doctest: +SKIP
    >>> print traj                                                 # doctest: +SKIP
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>  # doctest: +SKIP

    >>> traj2 = md.load_binpos('output.dcd', stride=2, top='topology.pdb') # doctest: +SKIP
    >>> print traj2                                                      # doctest: +SKIP
    <mdtraj.Trajectory with 250 frames, 423 atoms at 0x11136e410>         # doctest: +SKIP

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.BINPOSTrajectoryFile :  Low level interface to BINPOS files
    """
    from mdtraj.core.trajectory import _parse_topology

    # we make it not required in the signature, but required here. although this
    # is a little wierd, its good because this function is usually called by a
    # dispatch from load(), where top comes from **kwargs. So if its not supplied
    # we want to give the user an informative error message
    if top is None:
        raise ValueError('"top" argument is required for load_binpos')
    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_binpos. '
            'you supplied %s' % type(filename))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with BINPOSTrajectoryFile(str(filename)) as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(topology, n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)


cdef class BINPOSTrajectoryFile:
    """BINPOSTrajectoryFile(filename, mode='r', force_overwrite=True, **kwargs)

    Interface for reading and writing to an AMBER BINPOS file.
    This is a file-like object, that both reading or writing depending
    on the `mode` flag. It implements the context manager protocol,
    so you can also use it with the python 'with' statement.

    The conventional units in the BINPOS file are angstroms. The format only
    supports storing the cartesian coordinates.

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already
        exists on disk, should we overwrite it?

    Other Parameters
    ----------------
    min_chunk_size : int, default=100
        BINPOS, In read mode, we need to allocate a buffer in which to store the data
        without knowing how many frames are in the file. This parameter is the
        minimum size of the buffer to allocate.
    chunk_size_multiplier : int, default=1.5
        In read mode, we need to allocate a buffer in which to store the data
        without knowing how many frames are in the file. We can *guess* this
        information based on the size of the file on disk, but it's not perfect.
        This parameter inflates the guess by a multiplicative factor.

    Examples
    --------
    >>> # copy data from one file to another
    >>> with BINPOSTrajectoryFile('traj.binpos') as f:
    >>>     xyz = f.read()
    >>> with BINPOSTrajectoryFile('out.binpos') as f:
    >>>    f.write(xyz)

    See Also
    --------
    mdtraj.load_binpos : High-level wrapper that returns a ``md.Trajectory``
    """

    cdef int n_atoms
    cdef int n_frames # numnber of frames in the file, cached
    cdef int approx_n_frames
    cdef int min_chunk_size
    cdef int chunk_size_multiplier
    cdef int is_open
    cdef int force_overwrite
    cdef long int frame_counter
    cdef char* mode
    cdef char* filename
    cdef int write_initialized
    cdef readonly char* distance_unit

    cdef void* fh
    cdef molfile_timestep_t* timestep

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True, **kwargs):
        """Open an AMBER BINPOS file for reading/writing.
        """
        self.distance_unit = 'angstroms'
        self.is_open = False
        self.frame_counter = 0
        self.n_frames = -1
        self.force_overwrite = force_overwrite

        if str(mode) == 'r':
            self.n_atoms = 0
            if not os.path.exists(filename):
                raise IOError("The file '%s' doesn't exist" % filename)
            self.fh = open_binpos_read(filename, "binpos", &self.n_atoms)
            if self.n_atoms <= 0:
                raise IOError('Malformed BINPOS file. Number of atoms <= 0. '
                              'Are you sure this is a valid AMBER BINPOS file?')

            # binpos stores the 4 bytes per float, 3*n_atoms floats per frame, directly
            # with no compression
            self.approx_n_frames = os.stat(filename).st_size / (self.n_atoms * 4 * 3)

            self.min_chunk_size = max(kwargs.pop('min_chunk_size', 100), 1)
            self.chunk_size_multiplier = max(kwargs.pop('chunk_size_multiplier', 1.1), 0.01)
            self.is_open = True

            if self.fh is NULL:
                raise IOError('There was an error opening the binpos file: %s' % filename)

        elif str(mode) == 'w':
            self.write_initialized = False
        else:
            raise ValueError('mode must be one of "r" or "w". you supplied "%s"' % mode)

        self.timestep = <molfile_timestep_t*> malloc(sizeof(molfile_timestep_t))
        if self.timestep is NULL:
            raise MemoryError('There was an error allocating working space')

        for key in kwargs.keys():
            warnings.warn('kwarg "%s" was not recognized or processed' % key)

        self.filename = filename
        self.mode = mode

    def _initialize_write(self, n_atoms):
        force_overwrite = self.force_overwrite
        filename = self.filename
        if force_overwrite and os.path.exists(filename):
            os.unlink(filename)
        if not force_overwrite and os.path.exists(filename):
            raise IOError('"%s" already exists' % filename)
        self.fh = open_binpos_write(filename, "binpos", n_atoms)
        self.n_atoms = n_atoms
        self.is_open = True

    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """read_as_traj(topology, n_frames=None, stride=None, atom_indices=None)

        Read a trajectory from a BINPOS file

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

        initial = int(self.frame_counter)
        xyz = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        in_units_of(xyz, self.distance_unit, Trajectory._distance_unit, inplace=True)

        if stride is None:
            stride = 1
        time = (stride*np.arange(len(xyz))) + initial

        return Trajectory(xyz=xyz, topology=topology, time=time)

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """read(n_frames=None, stride=None, atom_indices=None)

        Read data from a BINPOS file

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : np.ndarray, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from
            the file. This may be slightly slower than the standard read
            because it required an extra copy, but will save memory.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
            The cartesian coordinates, in angstroms
        """
        if not str(self.mode) == 'r':
            raise ValueError('read() is only available when file is opened in mode="r"')

        if n_frames is not None:
            # if they supply the number of frames they want, that's easy
            if not int(n_frames) == n_frames:
                raise ValueError('n_frames must be an int, you supplied "%s"' % n_frames)
            return self._read(int(n_frames), atom_indices)[::stride]

        else:
            all_xyz = []
            while True:
                # guess the size of the chunk to read, based on how many frames we think are in the file
                # and how many we've currently read
                chunk = max(abs(int((self.approx_n_frames - self.frame_counter) * self.chunk_size_multiplier)),
                            self.min_chunk_size)
                xyz = self._read(chunk, atom_indices)

                if len(xyz) <= 0:
                    break
                all_xyz.append(xyz)

            return np.concatenate(all_xyz)[::stride]


    def _read(self, int n_frames, atom_indices=None):
        cdef int n_atoms_to_read, i

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

        # only used if atom_indices is given
        cdef np.ndarray[dtype=np.float32_t, ndim=2] framebuffer = np.zeros((self.n_atoms, 3), dtype=np.float32)

        cdef np.ndarray[dtype=np.float32_t, ndim=3] xyz = np.zeros((n_frames, n_atoms_to_read, 3), dtype=np.float32)

        for i in range(n_frames):
            if atom_indices is None:
                self.timestep.coords = &xyz[i,0,0]
                status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
            else:
                self.timestep.coords = &framebuffer[0, 0]
                status = read_next_timestep(self.fh, self.n_atoms, self.timestep)
                xyz[i, :, :] = framebuffer[atom_indices, :]

            self.frame_counter += 1

            if status != _BINPOS_SUCESS:
                break

        if status != _BINPOS_SUCESS:
            return xyz[:i]
        return xyz

    def write(self, xyz):
        """write(xyz)

        Write cartesian coordinates to a binpos file

        Parameters
        ----------
        xyz : np.ndarray, dtype=np.float32, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms in every frame, in angstroms.
        """
        xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False,
                          add_newaxis_on_deficient_ndim=True)

        if not self.write_initialized:
            self._initialize_write(xyz.shape[1])
            self.write_initialized = True
        else:
            if not self.n_atoms == xyz.shape[1]:
                raise ValueError('number of atoms in file (%d) does not match number '
                                 'of atoms youre trying to save (%d)' % (self.n_atoms, xyz.shape[1]))

        self._write(xyz)

    def _write(self, np.ndarray[dtype=np.float32_t, ndim=3] xyz not None):
        cdef int i, status
        cdef int n_frames = len(xyz)

        for i in range(n_frames):
            self.timestep.coords = &xyz[i, 0, 0]
            status = write_timestep(self.fh, self.timestep)

            if status != _BINPOS_SUCESS:
                raise RuntimeError("BINPOS Error: %s" % status)

    def seek(self, int offset, int whence=0):
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
        cdef int origin
        if str(self.mode) != 'r':
            raise NotImplementedError("seek is only supported in mode='r'")

        if whence == 0:
            origin = SEEK_SET
        elif whence == 1:
            origin = SEEK_CUR
        elif whence == 2:
            origin = SEEK_END
        else:
            raise IOError('Invalid argument')

        seek_timestep(self.fh, offset, origin)
        self.frame_counter = tell_timestep(self.fh)

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        return tell_timestep(self.fh)

    def close(self):
        """Close the BINPOS file"""
        if self.is_open:
            if str(self.mode) == 'r':
                close_file_read(self.fh)
            else:
                close_file_write(self.fh)
            self.is_open = False

    def __dealloc__(self):
        free(self.timestep)
        self.close()

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

    def __len__(self):
        if str(self.mode) != 'r':
            raise NotImplementedError('len() only available in mode="r" currently')
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.n_frames == -1:
            position = self.tell()
            self.seek(0, 2)
            self.n_frames = self.tell()
            self.seek(position)
        return self.n_frames
FormatRegistry.register_fileobject('.binpos')(BINPOSTrajectoryFile)
