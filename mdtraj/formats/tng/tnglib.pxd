#ctypedef long int64_t
from libc.stdint cimport int64_t

ctypedef enum tng_variable_n_atoms_flag: TNG_CONSTANT_N_ATOMS, TNG_VARIABLE_N_ATOMS
ctypedef enum tng_function_status: TNG_SUCCESS, TNG_FAILURE, TNG_CRITICAL

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

    # read frames in chunks?
    tng_function_status tng_util_pos_read_range(
                 const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **positions,
                 int64_t *stride_length)
