#include "fancy_index.hpp"
#include "stdio.h"
#include "stdlib.h"


/**
 * malloc and return a new range from 0 to n. The caller assumes ownership
 * of the memory, which must be `free`d.
 *
 **/
static int*
range(int n)
{
    int i;
    int* r = (int*) malloc(n * sizeof(int));
    if (r == NULL) { fprintf(stderr, "malloc failure in file '%s' in line %i\n", __FILE__, __LINE__); exit(1); }

    for (i = 0; i < n; i++) {
        r[i] = i;
    }
    return r;
}


/**
 * Numpy-like "fancy indexing" of a 2D array
 *
 * The equivalent python code for this function is
 *
 * >>> out = A[indx, indy]
 * 
 * where
 *   A.shape == (nx, ny),
 *   indx.shape == (nindx,)
 *   indy.shape == (nindy,)
 *
 *
 * The reason for writing this function in C is that from cython, we can't
 * call the numpy fancy indexing without acquiring the GIL, which means it
 * can't be used inside a prange construct
 **/
void
fancy_index2d(const float* A, int nx, int ny,
              const int* indx, int nindx, const int* indy, int nindy,
              float* out)
{
    int i, ii, j, jj;
    int* indx_ = (int*) indx;
    int* indy_ = (int*) indy;

    if (indx == NULL) { indx_ = range(nx); nindx = nx; }
    if (indy == NULL) { indy_ = range(ny); nindy = ny; }

    for (ii = 0; ii < nindx; ii++) {
        i = indx_[ii];
        for (jj = 0; jj < nindy; jj++) {
            j = indy_[jj];
            out[ii*nindy + jj] = A[i*ny + j];
        }
    }
    
    if (indx == NULL) { free(indx_); }
    if (indy == NULL) { free(indy_); }
}

