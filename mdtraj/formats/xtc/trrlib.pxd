cdef extern from "include/xdrfile.h":
    ctypedef struct XDRFILE:
        pass

    XDRFILE* xdrfile_open (char * path, char * mode)
    int xdrfile_close (XDRFILE * xfp)
    ctypedef float matrix[3][3]
    ctypedef float rvec[3]


cdef extern from "include/xdrfile_trr.h":
    int read_trr_natoms(char *fn, int *natoms)
    int read_trr_nframes(char* fn, unsigned long *nframes)

    # Read one frame of an open xtc file. If either of x, v, f, box are
    # NULL the arrays will be read from the file but not used.
    int read_trr(XDRFILE *xd, int natoms, int *step, float *t, float* lambd,
        matrix box, rvec* x, rvec* v, rvec* f)

    # Write a frame to xtc file
    int write_trr(XDRFILE *xd, int natoms, int step, float t, float lambd,
        matrix box, rvec* x, rvec* v, rvec* f)

cdef extern from "include/trr_header.h":
    ctypedef struct t_trnheader:
        int v_size
        int f_size

    # header
    int do_trnheader(XDRFILE *xd, int bRead, t_trnheader *sh)
