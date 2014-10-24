cdef extern from "include/molfile_plugin.h":
    ctypedef struct molfile_timestep_t:
      float *coords  # coordinates of all atoms, arranged xyzxyzxyz
      float *velocities  # space for velocities of all atoms; same layout
      float A, B, C, alpha, beta, gamma
      double physical_time

    ctypedef size_t molfile_ssize_t

    ctypedef struct molfile_timestep_metadata:
      unsigned int count                  # total # timesteps; -1 if unknown
      unsigned int avg_bytes_per_timestep # bytes per timestep
      int has_velocities

cdef extern from "include/dtrplugin.hxx":

    void* open_file_read(const char *path, const char *filetype, int *natoms)
    int read_timestep2(void *v, molfile_ssize_t n, molfile_timestep_t *ts)
    molfile_ssize_t read_times(void *v, molfile_ssize_t start, molfile_ssize_t count, double *times)
    void close_file_read(void *v)

    void* open_file_write(const char *path, const char *filetype, const int natoms)
    int write_timestep(void *v, const molfile_timestep_t *ts)
    void close_file_write(void *v)
    int read_timestep_metadata(void *v, molfile_timestep_metadata *m)
    ssize_t dtr_curframe(void* v)

