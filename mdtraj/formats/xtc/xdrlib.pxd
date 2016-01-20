cimport numpy as np
ctypedef np.npy_int64 int64_t

cdef extern from "include/xdrfile.h":
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    ctypedef float matrix[3][3]
    ctypedef float rvec[3]
    int xdrfile_close (XDRFILE * xfp)
    int xdr_seek(XDRFILE *xfp, int64_t pos, int whence)

cdef extern from "include/xdrfile_xtc.h":
    int read_xtc_natoms(char* fn, int* natoms)
    int read_xtc(XDRFILE *xd, int natoms, int *step, float *time, matrix box, rvec *x, float *prec)
    int write_xtc(XDRFILE *xd, int natoms, int step, float time, matrix box, rvec* x, float prec)
    int read_xtc_nframes(char* fn, unsigned long *nframes, unsigned long *est_nframes, int64_t** offsets)

