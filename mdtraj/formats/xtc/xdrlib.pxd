
cdef extern from "include/xdrfile.h":
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    ctypedef float matrix[3][3]
    ctypedef float rvec[3]
    int xdrfile_close (XDRFILE * xfp)

cdef extern from "include/xdrfile_xtc.h":
    int read_xtc_natoms(char* fn, int* natoms)
    int read_xtc(XDRFILE *xd, int natoms, int *step, float *time, matrix box, rvec *x, float *prec)
    int write_xtc(XDRFILE *xd, int natoms, int step, float time, matrix box, rvec* x, float prec)
    int xdrfile_read_int(int * ptr, int ndata, XDRFILE *xfp)


cimport numpy as np

ctypedef np.npy_int64 int64_t

cdef extern from "include/xdr_seek.h":
    int64_t xdr_tell(XDRFILE *xd)
    int xdr_seek(XDRFILE *xd, int64_t pos, int whence)
    int xdr_flush(XDRFILE* xd)
