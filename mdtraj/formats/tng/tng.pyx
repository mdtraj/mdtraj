# cython: c_string_type=str, c_string_encoding=ascii

cimport tnglib
from libc.stdio cimport printf
from libc.stdint cimport int64_t

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np
np.import_array()


cdef class TNGTrajectoryFile:
    cdef tnglib.tng_trajectory_t _traj
    cdef const char * filename
    cdef char mode
    cdef int is_open
    cdef readonly char * distance_unit
    cdef int64_t n_atoms  # number of atoms in the file
    cdef int64_t tot_n_frames
    cdef int64_t _pos

    def __cinit__(self, char * filename, char mode='r', force_overwrite=True, **kwargs):
        self.filename = filename
        self.mode = mode
        if mode == 'w':
            raise NotImplementedError()

        # TODO: TNG_CONSTANT_N_ATOMS assert that this is set

        if tnglib.tng_util_trajectory_open(filename, mode, & self._traj) == 0:
            self.is_open = True
        else:
            raise Exception("something went wrong during opening.")

        if tnglib.tng_num_particles_get(self._traj, & self.n_atoms) != 0:
            raise Exception("something went wrong during obtaining num particles")
        
        if tnglib.tng_num_frames_get(self._traj, & self.tot_n_frames) != 0:
            raise Exception("error during len determination")
        
        self._pos = 0

    def __len__(self):
        return self.tot_n_frames

    def __dealloc__(self):
        self.close()

    def close(self):
        "Close the XTC file handle"
        if self.is_open:
            tnglib.tng_util_trajectory_close(& self._traj)
            self.is_open = False
            self._pos = 0

    def _read(self, int n_frames, atom_indices):
        if self._pos >= self.tot_n_frames:
            raise EOFError()
        # buffer malloced by tng util function (tng_util_pos_read_range)
        cdef float * positions = NULL
        cdef int64_t stride_length, i, j
        # Set a default frame range
        cdef int64_t first_frame = self._pos, last_frame = self._pos + n_frames
        cdef int k
        
        if last_frame > self.tot_n_frames:# - 1):
            last_frame = self.tot_n_frames #- 1
        
        n_frames = last_frame - first_frame #+ 1
        
        #printf("_read: first_frame = %i; last_frame=%i\n", first_frame, last_frame)
        
        # handle atom_indices
        cdef int n_atoms_to_read
        if atom_indices is None:
            n_atoms_to_read = self.n_atoms
        elif isinstance(atom_indices, slice):
            n_atoms_to_read = len(np.arange(self.n_atoms)[atom_indices])

        cdef np.ndarray[ndim = 3, dtype = np.float32_t, mode = 'c'] xyz = \
            np.empty((n_frames, n_atoms_to_read, 3), dtype=np.float32)

        # Get the positions of all particles in the requested frame range.
        # The positions are stored in the positions array.
        # N.B. No proper error checks.
        # TODO: read time, step, box quantities
        if (tnglib.tng_util_pos_read_range(self._traj, self._pos, last_frame,
                                           & positions, & stride_length)
                == 0):
            # Print the positions of the wanted particle (zero based)
            for i in range(0, n_frames, stride_length):
                for j in range(0, self.n_atoms):
                    for k in range(3):
                        xyz[i, j, k] = positions[i / stride_length * self.n_atoms * 3
                                        + j * 3 + k]
        else:
            if positions: free(positions)
            raise Exception("Cannot read positions")
        if positions:
            free(positions)
            
        self._pos += n_frames

        return xyz

#     def read(self, n_frames=None, stride=None, atom_indices=None):
#         """read(n_frames=None, stride=None, atom_indices=None)
# 
#         Read data from a TNG file
# 
#         Parameters
#         ----------
#         n_frames : int, None
#             The number of frames you would like to read from the file.
#             If None, all of the remaining frames will be loaded.
#         stride : int, optional
#             Read only every stride-th frame.
#         atom_indices : array_like, optional
#             If not none, then read only a subset of the atoms coordinates from the
#             file. This may be slightly slower than the standard read because it required
#             an extra copy, but will save memory.
# 
#         Returns
#         -------
#         xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
#             The cartesian coordinates, in nanometers
#         time : np.ndarray, shape=(n_frames), dtype=np.float32
#             The simulation time, in picoseconds, corresponding to each frame
#         step : np.ndarray, shape=(n_frames), dtype=np.int32
#             The step in the simulation corresponding to each frame
#         box : np.ndarray, shape=(n_frames, 3, 3), dtype=np.float32
#             The box vectors in each frame.
# 
#         See Also
#         --------
#         read_as_traj : Returns a Trajectory object
#         """
#         if not str(self.mode) == 'r':
#             raise ValueError(
#                 'read() is only available when file is opened in mode="r"')
#         if not self.is_open:
#             raise IOError('file must be open to read from it.')
# 
#         if n_frames is not None:
#             # if they supply the number of frames they want, that's easy
#             if not int(n_frames) == n_frames:
#                 raise ValueError(
#                     'n_frames must be an int, you supplied "%s"' % n_frames)
#             xyz, time, step, box = self._read(int(n_frames), atom_indices)
#             xyz, time, step, box = xyz[::stride], time[
#                 ::stride], step[::stride], box[::stride]
#             if np.all(np.logical_and(box < 1e-10, box > -1e-10)):
#                 box = None
#             return xyz, time, step, box
# 
#         # if they want ALL of the remaining frames, we need to guess at the
#         # chunk size, and then check the exit status to make sure we're really
#         # at the EOF
#         all_xyz, all_time, all_step, all_box = [], [], [], []

#         while True:
#             # guess the size of the chunk to read, based on how many frames we
#             # think are in the file and how many we've currently read
#             chunk = max(abs(int((self.approx_n_frames - self.frame_counter) * self.chunk_size_multiplier)),
#                         self.min_chunk_size)
# 
#             xyz, time, step, box = self._read(chunk, atom_indices)
#             if len(xyz) <= 0:
#                 break
# 
#             all_xyz.append(xyz)
#             all_time.append(time)
#             all_step.append(step)
#             all_box.append(box)
# 
#         if len(all_xyz) == 0:
#             return np.array([]), np.array([]), np.array([]), np.array([])
#         all_xyz = np.concatenate(all_xyz)[::stride]
#         all_time = np.concatenate(all_time)[::stride]
#         all_step = np.concatenate(all_step)[::stride]
#         all_box = np.concatenate(all_box)[::stride]
#         if np.all(np.logical_and(all_box < 1e-10, all_box > -1e-10)):
#             all_box = None
#         return all_xyz, all_time, all_step, all_box
# 
#     def _read(self, int n_frames, atom_indices):
#         """Read a specified number of XTC frames from the buffer"""
# 
#         cdef int i = 0
#         cdef int status = _EXDROK
#         cdef int n_atoms_to_read
# 
#         if atom_indices is None:
#             n_atoms_to_read = self.n_atoms
#         elif isinstance(atom_indices, slice):
#             n_atoms_to_read = len(np.arange(self.n_atoms)[atom_indices])
#         else:
#             atom_indices = np.asarray(atom_indices)
#             if min(atom_indices) < 0:
#                 raise ValueError(
#                     'atom_indices should be zero indexed. you gave an index less than zero')
#             if max(atom_indices) >= self.n_atoms:
#                 raise ValueError(
#                     'atom indices should be zero indexed. you gave an index bigger than the number of atoms')
#             n_atoms_to_read = len(atom_indices)
# 
#         cdef np.ndarray[ndim = 3, dtype = np.float32_t, mode = 'c'] xyz = \
#             np.empty((n_frames, n_atoms_to_read, 3), dtype=np.float32)
#         cdef np.ndarray[ndim = 1, dtype = np.float32_t, mode = 'c'] time = \
#             np.empty((n_frames), dtype=np.float32)
#         cdef np.ndarray[ndim = 1, dtype = np.int32_t, mode = 'c'] step = \
#             np.empty((n_frames), dtype=np.int32)
#         cdef np.ndarray[ndim = 3, dtype = np.float32_t, mode = 'c'] box = \
#             np.empty((n_frames, 3, 3), dtype=np.float32)
#         cdef np.ndarray[ndim = 1, dtype = np.float32_t, mode = 'c'] prec = \
#             np.empty((n_frames), dtype=np.float32)
# 
#         # only used if atom_indices is given
#         cdef np.ndarray[dtype = np.float32_t, ndim = 2] framebuffer = np.zeros((self.n_atoms, 3), dtype=np.float32)
# 
#         while (i < n_frames) and (status != _EXDRENDOFFILE):
#             if atom_indices is None:
#                 status = xdrlib.read_xtc(self.fh, self.n_atoms, < int * > & step[i],
#                                          & time[i], < xdrlib.matrix > &box[i, 0, 0], < xdrlib.rvec * > & xyz[i, 0, 0], & prec[i])
#             else:
#                 status = xdrlib.read_xtc(self.fh, self.n_atoms, < int * > & step[i],
#                                          & time[i], < xdrlib.matrix > &box[i, 0, 0], < xdrlib.rvec * > & framebuffer[0, 0], & prec[i])
#                 xyz[i, :, :] = framebuffer[atom_indices, :]
# 
#             if status != _EXDRENDOFFILE and status != _EXDROK:
#                 raise RuntimeError('XTC read error: %s' %
#                                    _EXDR_ERROR_MESSAGES.get(status, 'unknown'))
#             i += 1
# 
#         if status == _EXDRENDOFFILE:
#             xyz = xyz[:i - 1]
#             box = box[:i - 1]
#             time = time[:i - 1]
#             step = step[:i - 1]
# 
#         self.frame_counter += i
# 
#         return xyz, time, step, box

    def __enter__(self):
        "Support the context manager protocol"
        return self

    def __exit__(self, *exc_info):
        "Support the context manager protocol"
        self.close()
