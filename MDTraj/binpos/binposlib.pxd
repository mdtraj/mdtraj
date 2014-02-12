cdef extern from "include/binposplugin.h":
    
    void *open_binpos_read(char *path, char *filetype, int *natoms)
    int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts)
    void close_file_read(void *v)

    void *open_binpos_write(char *path, char *filetype, int natoms)
    int write_timestep(void *v, molfile_timestep_t *ts)
    void close_file_write(void *v)
    int seek_timestep(void* v, long int offset, int origin)
    long int tell_timestep(void* v)

cdef extern from "include/molfile_plugin.h":
    ctypedef struct molfile_timestep_t:
      float *coords  # coordinates of all atoms, arranged xyzxyzxyz
      float *velocities  # space for velocities of all atoms; same layout
      float A, B, C, alpha, beta, gamma

