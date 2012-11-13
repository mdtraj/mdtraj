cdef extern from "include/xdrfile.h":
    ctypedef struct XDRFILE:
        pass


    XDRFILE* xdrfile_open (char * path, char * mode)
    int xdrfile_close (XDRFILE * xfp)
    int read_xtc_natoms(char* fn, int* natoms)
    int read_xtc(XDRFILE *xd, int natoms, int *step, float *time, float* box, float *x, float *prec)

    #int write_xtc(XDRFILE *xd, int natoms,int step,float time, matrix box, rvec *x, float prec)