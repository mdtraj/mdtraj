
cdef extern from "include/dcdplugin.h":
    ctypedef struct dcdhandle:
      long fd
      int natoms
      int nsets
      int setsread
      int istart
      int nsavc
      double delta
      int nfixed
      float *x, *y, *z
      int *freeind
      float *fixedcoords
      int reverse
      int charmm
      int first
      int with_unitcell

    dcdhandle* open_dcd_read(char *path, char *filetype, int *natoms, int* nsets)
    int read_next_timestep(dcdhandle *v, int natoms, molfile_timestep_t *ts)
    void close_file_read(dcdhandle *v)

    dcdhandle* open_dcd_write(const char *path, const char *filetype, const int natoms, const int with_unitcell)
    int write_timestep(dcdhandle *v, molfile_timestep_t *ts)
    void close_file_write(dcdhandle *v)
    int dcd_nsets(dcdhandle* v)
    int dcd_rewind(dcdhandle* dcd)


cdef extern from "include/molfile_plugin.h":
    ctypedef struct molfile_timestep_t:
      float *coords  # coordinates of all atoms, arranged xyzxyzxyz
      float *velocities  # space for velocities of all atoms; same layout
      float A, B, C, alpha, beta, gamma

