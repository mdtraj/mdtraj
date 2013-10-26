cdef extern from "include/xdrfile.h":
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    int xdrfile_close (XDRFILE * xfp)


cdef extern from "include/xdrfile_trr.h":

    int read_trr_natoms(char *fn, int *natoms)

    # Read one frame of an open xtc file. If either of x, v, f, box are
    # NULL the arrays will be read from the file but not used.
    int read_trr(XDRFILE *xd, int natoms, int *step, float *t, float* lambd,
        float* box, float* x, float* v, float* f)

    # Write a frame to xtc file
    int write_trr(XDRFILE *xd, int natoms, int step, float t, float lambd,
        float* box, float* x, float* v, float* f)
