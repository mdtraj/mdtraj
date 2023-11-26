# cython: c_string_type=str, c_string_encoding=ascii

from libc.stdio cimport printf
from libc.stdint cimport int64_t
from libc.stdlib cimport malloc, free

from math import ceil
from mdtraj.utils import ensure_type, cast_indices, in_units_of
from mdtraj.utils.six import string_types
from mdtraj.formats.registry import FormatRegistry
from mdtraj.formats import PDBTrajectoryFile
import mdtraj as md

import os
import numpy as np
cimport numpy as np
np.import_array()

__all__ = ['load_tng', 'TNGTrajectoryFile']

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
    tng_function_status tng_implicit_num_particles_set(
                const tng_trajectory_t tng_data,
                const int64_t n)

    # n frames
    tng_function_status tng_num_frames_get(
                const tng_trajectory_t tng_data,
                int64_t *n)
    
    # Units
    tng_function_status tng_distance_unit_exponential_get(
                const tng_trajectory_t tng_data,
                int64_t *exp);

    # Positions
    tng_function_status tng_util_pos_read_range(
                const tng_trajectory_t tng_data,
                const int64_t first_frame,
                const int64_t last_frame,
                float **positions,
                int64_t *stride_length)
    tng_function_status tng_util_pos_write(
                const tng_trajectory_t tng_data,
                const int64_t frame_nr,
                const float *positions)
    tng_function_status tng_util_pos_write_interval_set(
                const tng_trajectory_t tng_data,
                const int64_t i)

    # Box shape
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
    tng_function_status tng_util_box_shape_write(
                const tng_trajectory_t tng_data,
                const int64_t frame_nr,
                const float *box_shape)
    tng_function_status tng_util_box_shape_write_interval_set(
                const tng_trajectory_t tng_data,
                const int64_t i)

    # Time
    tng_function_status tng_util_time_of_frame_get(
                const tng_trajectory_t tng_data,
                const int64_t frame_nr,
                double *time)
    tng_function_status tng_time_per_frame_set(
                const tng_trajectory_t tng_data,
                const double time)
    tng_function_status tng_frame_set_first_frame_time_set(
                const tng_trajectory_t tng_data,
                const double first_frame_time)
    
    # Miscellaneous
    tng_function_status tng_num_frames_per_frame_set_get(
                const tng_trajectory_t tng_data,
                int64_t *n)
    tng_function_status tng_last_program_name_set(
                const tng_trajectory_t tng_data,
                const char *new_name)
    
    # Topology
    tng_function_status tng_chain_name_of_particle_nr_get(
                const tng_trajectory_t tng_data,
                const int64_t nr,
                char *name,
                const int max_len)
    tng_function_status tng_residue_name_of_particle_nr_get(
                const tng_trajectory_t tng_data,
                const int64_t nr,
                char *name,
                const int max_len)
    tng_function_status tng_global_residue_id_of_particle_nr_get(
                const tng_trajectory_t tng_data,
                const int64_t nr,
                int64_t *id)
    tng_function_status tng_atom_name_of_particle_nr_get(
                const tng_trajectory_t tng_data,
                const int64_t nr,
                char *name,
                const int max_len)
    tng_function_status tng_atom_type_of_particle_nr_get(
                const tng_trajectory_t tng_data,
                const int64_t nr,
                char *type,
                const int max_len)
    tng_function_status tng_molsystem_bonds_get(
                const tng_trajectory_t tng_data,
                int64_t *n_bonds,
                int64_t **from_atoms,
                int64_t **to_atoms)


@FormatRegistry.register_loader('.tng')
def load_tng(filename, top=None, stride=None, atom_indices=None, frame=None):
    """load_tng(filename, top=None, stride=None, atom_indices=None, frame=None)

    Load a Gromacs TNG file from disk.

    The .tng format is a cross-platform compressed binary trajectory format
    produced by the Gromacs software that stores atomic coordinates, box vectors,
    and time information. It optionally can also store topology information.
    It is lossy (storing coordinates to about 1e-3 nm) and extremely space-efficient.

    Parameters
    ----------
    filename : path-like
        Filename (string) of tng trajectory.
    top : {str, Trajectory, Topology}, optional
        If the TNG file does not contain topology information, it must be provided
        with this argument. Pass in either the path to a RCSB PDB file, a trajectory,
        or a topology to supply this information.
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

    if not isinstance(filename, (string_types, os.PathLike)):
        raise TypeError('filename must be of type path-like for load_tng. '
                        'you supplied %s' % type(filename))
    if top is not None:
        topology = _parse_topology(top)
    else:
        topology = None
    atom_indices = cast_indices(atom_indices)

    with TNGTrajectoryFile(str(filename), 'r') as f:
        if frame is not None:
            f.seek(frame)
            n_frames = 1
        else:
            n_frames = None

        return f.read_as_traj(topology, n_frames=n_frames, stride=stride,
                              atom_indices=atom_indices)

cdef class TNGTrajectoryFile(object):
    """TNGTrajectoryFile(filename, mode='r', force_overwrite=True, **kwargs)

    Interface for reading and writing to a Gromacs TNG file.
    This is a file-like object that supports both reading and writing.
    It also supports the context manager protocol, so you can use it
    with the python 'with' statement.

    This class supports both reading and writing of coordinates, time,
    and unit cell parameters.  It also supports reading (but not writing)
    of topology information.

    Parameters
    ----------
    filename : str
        The filename to open. A path to a file on disk.
    mode : {'r', 'w'}
        The mode in which to open the file, either 'r' for read or 'w' for write.
    force_overwrite : bool
        If opened in write mode, and a file by the name of `filename` already exists on disk, should we overwrite it?

    Examples
    --------
    >>> # read the data from from an TNG file
    >>> with TNGTrajectoryFile('traj.xtc') as f:
    >>>    xyz, time, box = f.read()
    >>>    top = f.topology

    >>> # write some random coordinates to an TNG file
    >>> with TNGTrajectoryFile('output.xtc', 'w') as f:
    >>>     f.write(np.random.randn(10,1,3))

    See Also
    --------
    mdtraj.load_tng : High-level wrapper that returns a ``md.Trajectory``
    """
    cdef tng_trajectory_t _traj
    cdef const char * filename
    cdef char mode
    cdef int is_open
    cdef float _distance_scale
    cdef int64_t n_atoms  # number of atoms in the file
    cdef int64_t tot_n_frames
    cdef int64_t _pos
    cdef int64_t _frames_per_frame_set
    cdef float _time_per_frame
    cdef object _topology
    cdef readonly char* distance_unit

    def __cinit__(self, char* filename, char* mode='r', force_overwrite=True, **kwargs):
        self.distance_unit = 'nanometers'
        self.filename = filename
        self.mode = mode[0]
        if self.mode == b'w' and os.path.exists(filename):
            if force_overwrite:
                os.unlink(filename)
            else:
                raise IOError('"%s" already exists' % filename)
        if self.mode == b'r' and not os.path.exists(filename):
            raise IOError('"%s" does not exist' % filename)
        if tng_util_trajectory_open(filename, self.mode, & self._traj) == TNG_SUCCESS:
            self.is_open = True
        else:
            raise Exception("An error ocurred opening the file.")
        cdef int64_t exponent;
        if self.mode == b'r':
            if tng_num_particles_get(self._traj, & self.n_atoms) != TNG_SUCCESS:
                raise Exception("something went wrong during obtaining num particles")
            
            if tng_num_frames_get(self._traj, & self.tot_n_frames) != TNG_SUCCESS:
                raise Exception("error during len determination")
     
            if tng_distance_unit_exponential_get(self._traj, &exponent) != TNG_SUCCESS:
                raise Exception("Error reading distance unit exponent")
            self._distance_scale = 10.0**(exponent+9)
        else:
            self.tot_n_frames = 0
            self._distance_scale = 1.0
            if tng_num_frames_per_frame_set_get(self._traj, &self._frames_per_frame_set) != TNG_SUCCESS:
                raise Exception("Error reading number of frames per frame set")
            if tng_last_program_name_set(self._traj, 'MDTraj %s' % md.version.version) != TNG_SUCCESS:
                raise Exception("Error writing program name")
        self._pos = 0
        self._time_per_frame = 0
        self._topology = self._read_topology()

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
    
    def _read_topology(self):
        cdef char text[1024]
        cdef int64_t residue_id
        last_chain_name = None
        last_residue_id = None
        top = md.Topology()
        PDBTrajectoryFile._loadNameReplacementTables()
        
        # Loop over atoms, load the information for each one, and create the Topology
        
        cdef int i
        for i in range(self.n_atoms):
            if tng_chain_name_of_particle_nr_get(self._traj, i, text, 1024) != TNG_SUCCESS:
                return None
            chain_name = str(text)
            if tng_residue_name_of_particle_nr_get(self._traj, i, text, 1024) != TNG_SUCCESS:
                return None
            residue_name = str(text).strip()
            if residue_name in PDBTrajectoryFile._residueNameReplacements:
                residue_name = PDBTrajectoryFile._residueNameReplacements[residue_name]
            if tng_global_residue_id_of_particle_nr_get(self._traj, i, &residue_id) != TNG_SUCCESS:
                return None
            if tng_atom_name_of_particle_nr_get(self._traj, i, text, 1024) != TNG_SUCCESS:
                return None
            atom_name = str(text).strip()
            if residue_name in PDBTrajectoryFile._atomNameReplacements and atom_name in PDBTrajectoryFile._atomNameReplacements[residue_name]:
                atom_name = PDBTrajectoryFile._atomNameReplacements[residue_name][atom_name]
            if tng_atom_type_of_particle_nr_get(self._traj, i, text, 1024) != TNG_SUCCESS:
                return None
            try:
                element = md.element.get_by_symbol(str(text))
            except KeyError:
                element = None
            if chain_name != last_chain_name:
                chain = top.add_chain()
                last_chain_name = chain_name
                last_residue_id = None
            if residue_id != last_residue_id:
                residue = top.add_residue(residue_name, chain, residue_id)
                last_residue_id = residue_id
            top.add_atom(atom_name, element, residue)

        if all(r.name == '' for r in top.residues) and all(a.name == '' for a in top.atoms):
            # This file doesn't contain topology information.
            return None

        # Now that we know how many atoms are in each residue, we can guess the elements.

        for atom in top.atoms:
            if atom.element is None:
                atom.element = PDBTrajectoryFile._guess_element(atom.name, atom.residue.name, atom.residue.n_atoms)
        
        # Load bonds.
        
        cdef int64_t n_bonds
        cdef int64_t* from_atoms
        cdef int64_t* to_atoms
        if tng_molsystem_bonds_get(self._traj, &n_bonds, &from_atoms, &to_atoms) != TNG_SUCCESS:
            raise Exception("Error reading bonds")
        for i in  range(n_bonds):
            top.add_bond(top.atom(from_atoms[i]), top.atom(to_atoms[i]))
        free(from_atoms)
        free(to_atoms)
        return top

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
                xyz *= self._distance_scale
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
            time = frame_time*1e12

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
                box *= self._distance_scale
            else:
                raise Exception("Error reading box shape")
        finally:
            if box_shape != NULL:
                free(box_shape)

        self._pos += 1

        return xyz, time, box
    
    def read_as_traj(self, topology=None, n_frames=None, stride=None, atom_indices=None):
        """read_as_traj(topology, n_frames=None, stride=None, atom_indices=None)

        Read a trajectory from a TNG file

        Parameters
        ----------
        topology : Topology, default=None, optional
            The system topology. If None, the Topology read from the file will be used.
        n_frames : int, default=None
            The number of frames you would like to read from the file.
            If None, all of the remaining frames will be loaded.
        stride : int, default=None, optional
            Read only every stride-th frame.
        atom_indices : array_like, default=None, optional
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
        if topology is None:
            topology = self.topology
        if topology is None:
            raise ValueError('This file does not contain topology information, and no Topology was provided.')
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
        if not self.mode == b'r':
            raise ValueError(
                'read() is only available when file is opened in mode="r"')
        if not self.is_open:
            raise IOError('file must be open to read from it.')

        if stride is None:
            stride = 1
        max_frames = int(ceil((self.tot_n_frames-self._pos)/float(stride)))
        if n_frames is None or n_frames > max_frames:
            n_frames = max_frames
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

        all_xyz = np.empty((n_frames, len(atom_indices), 3), dtype=np.float32)
        all_time = []
        all_box = np.empty((n_frames, 3, 3), dtype=np.float32);
        for i in range(n_frames):
            xyz, time, box = self._read_frame(atom_indices)
            all_xyz[i] = xyz
            all_time.append(time)
            all_box[i] = box
            self._pos += stride-1
        if any(time is None for time in all_time):
            all_time = None
        else:
            all_time = np.array(all_time, dtype=np.float32)
        if np.all(np.logical_and(all_box < 1e-10, all_box > -1e-10)):
            all_box = None
        return all_xyz, all_time, all_box

    def write(self, xyz, time=None, box=None):
        """write(xyz, time=None, box=None)

        Write data to a TNG file

        Parameters
        ----------
        xyz : np.ndarray, dtype=np.float32, shape=(n_frames, n_atoms, 3)
            The cartesian coordinates of the atoms, in nanometers
        time : np.ndarray, dtype=float32, shape=(n_frames), optional
            The simulation time corresponding to each frame, in picoseconds.
        box : np.ndarray, dtype=float32, shape=(n_frames, 3, 3), optional
            The periodic box vectors of the simulation in each frame, in nanometers.
            If not supplied, all box vectors will be recorded as (0,0,0).
        """
        if self.mode != b'w':
            raise ValueError('write() is only available when the file is opened in mode="w"')

        # do typechecking, and then dispatch to the c level function
        xyz = ensure_type(xyz, dtype=np.float32, ndim=3, name='xyz', can_be_none=False,
                          add_newaxis_on_deficient_ndim=True, warn_on_cast=False)
        n_frames = len(xyz)
        time = ensure_type(time, dtype=np.float32, ndim=1, name='time', can_be_none=True,
                           shape=(n_frames,), add_newaxis_on_deficient_ndim=True,
                           warn_on_cast=False)
        box = ensure_type(box, dtype=np.float32, ndim=3, name='box', can_be_none=True,
                          shape=(n_frames, 3, 3), add_newaxis_on_deficient_ndim=True,
                          warn_on_cast=False)

        if self.tot_n_frames == 0:
            # Initialization before writing the first frame.
            self.n_atoms = xyz.shape[1]
            if tng_implicit_num_particles_set(self._traj, self.n_atoms) != TNG_SUCCESS:
                raise Exception("Error writing number of particles")
            if tng_util_pos_write_interval_set(self._traj, 1) != TNG_SUCCESS:
                raise Exception("Error setting position stride")
            if tng_util_box_shape_write_interval_set(self._traj, 1) != TNG_SUCCESS:
                raise Exception("Error setting box shape stride")
        else:
            if not self.n_atoms == xyz.shape[1]:
                raise ValueError("This file has %d atoms, but you're now "
                    "trying to write %d atoms" % (self.n_atoms, xyz.shape[1]))
        if n_frames > 1 and time is not None:
            # Compute the time per frame and make sure it's consistent.
            time_per_frame = time[1]-time[0]
            if (self._time_per_frame != 0 and self._time_per_frame != time_per_frame) or not np.allclose(time[:-1]+time_per_frame, time[1:]):
                raise ValueError("Frames must be evenly spaced in time")
            if self._time_per_frame == 0:
                self._time_per_frame = time_per_frame
                if tng_time_per_frame_set(self._traj, time_per_frame*1e-12) != TNG_SUCCESS:
                    raise Exception("Error writing time per frame")

        if box is None:
            # Make each box[i] be the all zeros, which indicates the lack of a unitcell
            box = np.zeros((n_frames, 3, 3), dtype=np.float32)

        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] xyz_array = xyz
        cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] box_array = box
        for i in range(len(xyz)):
            if tng_util_pos_write(self._traj, self.tot_n_frames, &xyz_array[i,0,0]) != TNG_SUCCESS:
                raise Exception("Error writing positions")
            if tng_util_box_shape_write(self._traj, self.tot_n_frames, &box_array[i,0,0]) != TNG_SUCCESS:
                raise Exception("Error writing box shape")
            if self.tot_n_frames%self._frames_per_frame_set == 0 and time is not None:
                if tng_frame_set_first_frame_time_set(self._traj, time[i]*1e-12) != TNG_SUCCESS:
                    raise Exception("Error writing first frame time")
            self.tot_n_frames += 1

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

        if self.mode != b'r':
            raise NotImplementedError('seek() only available in mode="r" currently')
        if whence == 0 and offset >= 0:
            self._pos = offset
        elif whence == 1:
            self._pos += offset
        elif whence == 2 and offset <= 0:
            self._pos = self.tot_n_frames+offset
        else:
            raise IOError('Invalid argument')

    def tell(self):
        """Current file position

        Returns
        -------
        offset : int
            The current frame in the file.
        """
        if self.mode != b'r':
            raise NotImplementedError('tell() only available in mode="r" currently')
        return int(self._pos)

    @property
    def topology(self):
        """The Topology loaded from this TNG file, or None if it does not contain topology information. Only available when the file is opened in mode='r'.
        """
        return self._topology

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()

FormatRegistry.register_fileobject('.tng')(TNGTrajectoryFile)
