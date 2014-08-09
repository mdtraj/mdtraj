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
    int read_xtc_nframes(char* fn, unsigned long *nframes)