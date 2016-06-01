# cython: c_string_type=str, c_string_encoding=ascii

from libc.stdio cimport printf
from libc.stdint cimport int64_t
from libc.stdlib cimport malloc, free

from math import ceil
from mdtraj.utils import cast_indices, in_units_of
from mdtraj.utils.six import string_types
from mdtraj.formats.registry import FormatRegistry

import numpy as np
cimport numpy as np
np.import_array()

ctypedef enum tng_variable_n_atoms_flag: TNG_CONSTANT_N_ATOMS, TNG_VARIABLE_N_ATOMS
ctypedef enum tng_function_status: TNG_SUCCESS, TNG_FAILURE, TNG_CRITICAL
ctypedef enum tng_data_type: TNG_CHAR_DATA, TNG_INT_DATA, TNG_FLOAT_DATA, TNG_DOUBLE_DATA
ctypedef enum tng_hash_mode: TNG_SKIP_HASH, TNG_USE_HASH

cdef long long TNG_TRAJ_BOX_SHAPE = 0x0000000010000000LL

cdef extern from "tng/tng_io.h":
    ctypedef struct tng_trajectory_t:
        pass

    # open/close
    tng_function_status tng_util_trajectory_open(
                const char *filename,
                const char mode,
                tng_trajectory_t *tng_data_p)
    tng_function_status tng_util_trajectory_close(
                tng_trajectory_t *tng_data_p)

    # n particles
    tng_function_status tng_num_particles_get(
                const tng_trajectory_t tng_data,
                int64_t *n)

    # n frames
    tng_function_status tng_num_frames_get(
                const tng_trajectory_t tng_data,
                int64_t *n)
    
    # Units
    tng_function_status tng_distance_unit_exponential_get(
                const tng_trajectory_t tng_data,
                int64_t *exp);

    # read frames in chunks
    tng_function_status tng_util_pos_read_range(
                const tng_trajectory_t tng_data,
                const int64_t first_frame,
                const int64_t last_frame,
                float **positions,
                int64_t *stride_length)

    # Used to read the box shape
    tng_function_status tng_data_vector_interval_get(
                const tng_trajectory_t tng_data,
                const int64_t block_id,
                const int64_t start_frame_nr,
                const int64_t end_frame_nr,
                const char hash_mode,
                void **values,
                int64_t *stride_length,
                int64_t *n_values_per_frame,
                char *type)

    # Read time
    tng_function_status tng_util_time_of_frame_get(
                const tng_trajectory_t tng_data,
                const int64_t frame_nr,
                double *time)


@FormatRegistry.register_loader('.tng')
def load_tng(filename, top=None, stride=None, atom_indices=None, frame=None):
    """load_tng(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load a Gromacs TNG file from disk.

    The .tng format is a cross-platform compressed binary trajectory format
    produced by the gromacs software that stores atomic coordinates, box
    vectors, and time information. It is lossy (storing coordinates to about
    1e-3 A) and extremely space-efficient.

    Parameters
    ----------
    filename : str
        Filename (string) of tng trajectory.
    top : {str, Trajectory, Topology}
        The TNG format does not contain topology information. Pass in either the
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
    >>> traj = md.load_tng('output.tng', top='topology.pdb')
    >>> print traj
    <mdtraj.Trajectory with 500 frames, 423 atoms at 0x110740a90>

    Returns
    -------
    trajectory : md.Trajectory
        The resulting trajectory, as an md.Trajectory object.

    See Also
    --------
    mdtraj.TNGTrajectoryFile :  Low level interface to TNG files
    """
    from mdtraj.core.trajectory import _parse_topology
    if top is None:
        raise ValueError('"top" argument is required for load_tng')

    if not isinstance(filename, string_types):
        raise TypeError('filename must be of type string for load_tng. '
                        'you supplied %s' % type(filename))

    topology = _parse_topology(top)
    atom_indices = cast_indices(atom_indices)

    with TNGTrajectoryFile(filename, 'r') as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(topology, n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)

cdef class TNGTrajectoryFile(object):
    cdef tng_trajectory_t _traj
    cdef const char * filename
    cdef char mode
    cdef int is_open
    cdef float distance_scale
    cdef int64_t n_atoms  # number of atoms in the file
    cdef int64_t tot_n_frames
    cdef int64_t _pos

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True, **kwargs):
        self.filename = filename
        self.mode = mode[0]
        if self.mode == 'w':
            raise NotImplementedError()

        # TODO: TNG_CONSTANT_N_ATOMS assert that this is set

        if tng_util_trajectory_open(filename, self.mode, & self._traj) == TNG_SUCCESS:
            self.is_open = True
        else:
            raise Exception("something went wrong during opening.")

        if tng_num_particles_get(self._traj, & self.n_atoms) != TNG_SUCCESS:
            raise Exception("something went wrong during obtaining num particles")
        
        if tng_num_frames_get(self._traj, & self.tot_n_frames) != TNG_SUCCESS:
            raise Exception("error during len determination")
 
        cdef int64_t exponent;
        if tng_distance_unit_exponential_get(self._traj, &exponent) != TNG_SUCCESS:
            raise Exception("Error reading distance unit exponent")
        self.distance_scale = 10.0**(exponent+9)
        self._pos = 0

    def __len__(self):
        return self.tot_n_frames

    def __dealloc__(self):
        self.close()

    def close(self):
        "Close the TNG file handle"
        if self.is_open:
            tng_util_trajectory_close(& self._traj)
            self.is_open = False
            self._pos = 0

    def _read_frame(self, atom_indices):
        if self._pos >= self.tot_n_frames:
            raise EOFError()
        cdef int64_t stride_length, i, j, n_values_per_frame
        # Set a default frame range
        cdef int k
        cdef char data_type
        
        # handle atom_indices
        cdef int n_atoms_to_read = len(atom_indices)

        # Get the positions of all particles in the requested frame range.
        
        cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] xyz = np.empty((n_atoms_to_read, 3), dtype=np.float32)
        cdef float* positions = NULL
        try:
            if tng_util_pos_read_range(self._traj, self._pos, self._pos, &positions, &stride_length) == TNG_SUCCESS:
                for i, atom_index in enumerate(atom_indices):
                    for j in range(3):
                        xyz[i, j] = positions[atom_index*3 + j]
                xyz *= self.distance_scale
            else:
                raise Exception("Error reading positions")
        finally:
            if positions != NULL:
                free(positions)

        # Get the time of each frame in the requested range.
        
        cdef double frame_time
        if tng_util_time_of_frame_get(self._traj, self._pos, &frame_time) != TNG_SUCCESS:
            # This file doesn't specify times.
            time = None
        else:
            time = frame_time

        # Get the box shape of each frame in the requested range.
        
        cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] box = np.empty((3, 3), dtype=np.float32)
        cdef void* box_shape = NULL
        cdef float* float_box
        cdef double* double_box
        try:
            if tng_data_vector_interval_get(self._traj, TNG_TRAJ_BOX_SHAPE, self._pos, self._pos, TNG_USE_HASH,
                                        &box_shape, &stride_length, &n_values_per_frame, &data_type) == TNG_SUCCESS:
                if data_type == TNG_DOUBLE_DATA:
                    double_box = <double*>box_shape
                    for j in range(3):
                        for k in range(3):
                            box[j, k] = double_box[j*3 + k]
                else:
                    float_box = <float*>box_shape
                    for j in range(3):
                        for k in range(3):
                            box[j, k] = float_box[j*3 + k]
                box *= self.distance_scale
            else:
                raise Exception("Error reading box shape")
        finally:
            if box_shape != NULL:
                free(box_shape)

        self._pos += 1

        return xyz, time, box
    def read_as_traj(self, topology, n_frames=None, stride=None, atom_indices=None):
        """read_as_traj(topology, n_frames=None, stride=None, atom_indices=None)

        Read a trajectory from a TNG file

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

        xyz, time, box = self.read(n_frames=n_frames, stride=stride, atom_indices=atom_indices)
        if len(xyz) == 0:
            return Trajectory(xyz=np.zeros((0, topology.n_atoms, 3)), topology=topology)

        trajectory = Trajectory(xyz=xyz, topology=topology, time=time)
        trajectory.unitcell_vectors = box
        return trajectory

    def read(self, n_frames=None, stride=None, atom_indices=None):
        """read(n_frames=None, stride=None, atom_indices=None)

        Read data from a TNG file

        Parameters
        ----------
        n_frames : int, None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : int, optional
            Read only every stride-th frame.
        atom_indices : array_like, optional
            If not none, then read only a subset of the atoms coordinates from the
            file. This may be slightly slower than the standard read because it requires
            an extra copy, but will save memory.

        Returns
        -------
        xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
            The cartesian coordinates, in nanometers
        time : np.ndarray, shape=(n_frames), dtype=np.float32
            The simulation time, in picoseconds, corresponding to each frame
        box : np.ndarray, shape=(n_frames, 3, 3), dtype=np.float32
            The box vectors in each frame.

        See Also
        --------
        read_as_traj : Returns a Trajectory object
        """
        if not self.mode == 'r':
            raise ValueError(
                'read() is only available when file is opened in mode="r"')
        if not self.is_open:
            raise IOError('file must be open to read from it.')

        if stride is None:
            stride = 1
        if n_frames is None:
            n_frames = ceil((self.tot_n_frames-self._pos)/float(stride))
        if atom_indices is None:
            atom_indices = np.arange(self.n_atoms)
        elif isinstance(atom_indices, slice):
            atom_indices = np.arange(self.n_atoms)[atom_indices]
        else:
            atom_indices = np.asarray(atom_indices)
            if min(atom_indices) < 0:
                raise ValueError('atom_indices should be zero indexed. You gave an index less than zero')
            if max(atom_indices) >= self.n_atoms:
                raise ValueError('atom indices should be zero indexed. You gave an index bigger than the number of atoms')
        
        # Read the frames one at a time.

        all_xyz = np.empty((n_frames, len(atom_indices), 3))
        all_time = []
        all_box = np.empty((n_frames, 3, 3));
        for i in range(n_frames):
            xyz, time, box = self._read_frame(atom_indices)
            all_xyz[i] = xyz
            all_time.append(time)
            all_box[i] = box
            self._pos += stride-1
        if any(time is None for time in all_time):
            all_time = None
        else:
            all_time = np.array(all_time)
        if np.all(np.logical_and(all_box < 1e-10, all_box > -1e-10)):
            all_box = None
        return all_xyz, all_time, all_box

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

        if self.mode != 'r':
            raise NotImplementedError('seek() only available in mode="r" currently')
        if whence == 0 and offset >= 0:
            self._pos = offset
        elif whence == 1:
            self._pos += offset
        elif whence == 2 and offset <= 0:
            self._pos = self._tot_n_frames-offset-1
        else:
            raise IOError('Invalid argument')

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()
