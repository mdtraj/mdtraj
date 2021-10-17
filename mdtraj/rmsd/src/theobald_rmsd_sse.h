#include "msvccompat.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "util_sse.h"
#include "theobald_rmsd.h"
#include "stdio.h"


float msd_axis_major(const int nrealatoms, const int npaddedatoms, const int rowstride,
                     const float* aT, const float* bT, const float G_a, const float G_b) {
    /*   Computes the mean-square-deviation between two centered structures in axis-major format.
     *
     *   Structure setup for this function:
     *
     *   structures are stored axis major, and if this file is compiled with 
     *   -DALIGNED, possibly with extra padding to ensure you meet two constraints:
     *       - the number of elements in a row must be a multiple of 4
     *       - the first element in each row must be aligned to a 16 byte boundary
     *
     *   note that if you meet the second condition for the first row, and meet the
     *   first condition, the alignment will automatically be satisfied for every row.
     *
     *   the layout in memory for a structure of 7 atoms would look like this (read row-major):
     *
     *       x0 x1 x2 x3 x4 x5 x6 0
     *       y0 y1 y2 y3 y4 y5 y6 0
     *       z0 z1 z2 z3 z4 z5 z6 0
     *
     *   if your structure has a number of atoms that is not a multiple of 4, you must
     *   pad it out to a multiple of 4 using zeros (using anything other than zero will
     *   make the calculation go wrong).
     *
     * On the other hand, when this file is compiled without -DALIGNED, then 
     * there are no 16 byte alignment or dummy atom requirements, and the
     * "npaddedatoms" argument is ignored.
     *
     *   arguments:
     *       nrealatoms:   the *actual* number of atoms in the structure
     *
     *       npaddedatoms: the number of atoms in the structure including padding atoms;
     *                     should equal nrealatoms rounded up to the next multiple of 4
     *                     THIS IS ONLY USED IF COMPILED WITH -DALIGNED. OTHERWISE
     *                     IT IS IGNORED.
     *
     *       rowstride:    the offset in elements between rows in the arrays. will prob
     *                     be equal to npaddedatoms, but you might use something else if
     *                     (for example) you were subsetting the structure
     *
     *       aT:           pointer to start of first structure (A). should be aligned to
     *                     a 16-byte boundary if compiled with -DALIGNED
     *
     *       bT:           pointer to start of second structure (B). should be aligned to
     *                     a 16-byte boundary if compiled with -DALIGNED
     *
     *       G_a:          trace of A'A
     *
     *       G_b:          trace of B'B
     */
    int niters, k;
#ifndef ALIGNED
    static const int masks[4][4] = {
        {1, 1, 1, 1},
        {1, 0, 0, 0},
        {1, 1, 0, 0},
        {1, 1, 1, 0}
    };
    int const *mask;
#endif
    __m128 xx,xy,xz,yx,yy,yz,zx,zy,zz;
    __m128 ax,ay,az,b;
    __m128 t0,t1,t2;
    /* Will have 3 garbage elements at the end */
    _ALIGNED(16) float M[12];
    const float* aTx = aT;
    const float* aTy = aT+rowstride;
    const float* aTz = aT+2*rowstride;
    const float* bTx = bT;
    const float* bTy = bT+rowstride;
    const float* bTz = bT+2*rowstride;

    if (aT==bT && G_a==G_b)
        return 0.0;

#ifdef ALIGNED
    niters = npaddedatoms >> 2;
    /* npaddedatoms must be a multiple of 4 */
    assert(npaddedatoms % 4 == 0);
#else
    niters = (nrealatoms + 4-1) / 4;
    mask = masks[nrealatoms%4];
#endif
    
    xx = xy = xz = yx = yy = yz = zx = zy = zz = _mm_setzero_ps();
    for (k = 0; k < niters; k++) {
#ifdef ALIGNED
        ax = _mm_load_ps(aTx);
        ay = _mm_load_ps(aTy);
        az = _mm_load_ps(aTz);
        b = _mm_load_ps(bTx);
#else
        if (k == niters - 1) {
            ax = _mm_set_ps(mask[0] ? aTx[0] : 0, mask[1] ? aTx[1] : 0, mask[2] ? aTx[2] : 0, mask[3] ? aTx[3] : 0);
            ay = _mm_set_ps(mask[0] ? aTy[0] : 0, mask[1] ? aTy[1] : 0, mask[2] ? aTy[2] : 0, mask[3] ? aTy[3] : 0);
            az = _mm_set_ps(mask[0] ? aTz[0] : 0, mask[1] ? aTz[1] : 0, mask[2] ? aTz[2] : 0, mask[3] ? aTz[3] : 0);
            b =  _mm_set_ps(mask[0] ? bTx[0] : 0, mask[1] ? bTx[1] : 0, mask[2] ? bTx[2] : 0, mask[3] ? bTx[3] : 0);
        } else {
            ax = _mm_loadu_ps(aTx);
            ay = _mm_loadu_ps(aTy);
            az = _mm_loadu_ps(aTz);
            b = _mm_loadu_ps(bTx);
        }
#endif

        t0 = ax;
        t1 = ay;
        t2 = az;

        t0 = _mm_mul_ps(t0,b);
        t1 = _mm_mul_ps(t1,b);
        t2 = _mm_mul_ps(t2,b);

        xx = _mm_add_ps(xx,t0);
        yx = _mm_add_ps(yx,t1);
        zx = _mm_add_ps(zx,t2);

#ifdef ALIGNED
        b = _mm_load_ps(bTy);
#else
        if (k == niters - 1) {
            b =  _mm_set_ps(mask[0] ? bTy[0] : 0, mask[1] ? bTy[1] : 0, mask[2] ? bTy[2] : 0, mask[3] ? bTy[3] : 0);
        } else {
            b = _mm_loadu_ps(bTy);
        }
#endif
        t0 = ax;
        t1 = ay;
        t2 = az;

        t0 = _mm_mul_ps(t0,b);
        t1 = _mm_mul_ps(t1,b);
        t2 = _mm_mul_ps(t2,b);

        xy = _mm_add_ps(xy,t0);
        yy = _mm_add_ps(yy,t1);
        zy = _mm_add_ps(zy,t2);

#ifdef ALIGNED
        b = _mm_load_ps(bTz);
#else
        if (k == niters - 1) {
            b =  _mm_set_ps(mask[0] ? bTz[0] : 0, mask[1] ? bTz[1] : 0, mask[2] ? bTz[2] : 0, mask[3] ? bTz[3] : 0);
        } else {
            b = _mm_loadu_ps(bTz);
        }
#endif

        ax = _mm_mul_ps(ax,b);
        ay = _mm_mul_ps(ay,b);
        az = _mm_mul_ps(az,b);

        xz = _mm_add_ps(xz,ax);
        yz = _mm_add_ps(yz,ay);
        zz = _mm_add_ps(zz,az);

        aTx += 4;
        aTy += 4;
        aTz += 4;
        bTx += 4;
        bTy += 4;
        bTz += 4;
    }
    /* Epilogue - reduce 4 wide vectors to one wide */
    REDUCTION_EPILOGUE(xx, xy, xz, yx, yy, yz, zx, zy, zz, t0, t1, t2);

    _mm_store_ps(M  , xx);
    _mm_store_ps(M+4, yy);
    _mm_store_ps(M+8, zz);

    return msdFromMandG(M, G_a, G_b, nrealatoms, 0, NULL);
}


float msd_atom_major(const int nrealatoms, const int npaddedatoms,
                     const float* a, const float* b, const float G_a, const float G_b,
                     int computeRot, float rot[9]) {
    /* Computes the mean-square-deviation between two centered structures in
     * atom-major format.
     *
     * Structure setup for this function:
     *
     *   If this file is compiled with -DALIGNED, structures are stored atom
     *   major obeying two constraints:
     *       - if the number of atoms is not divisible by four, the structure is padded out
     *         with dummy atoms with zero in each coordinate up to an even multiple of 4 atoms.
     *       - the structure is aligned to a 16-byte boundary
     *
     *   the layout in memory for a structure of 7 atoms would look like this (read row-major):
     *
     *       x0 y0 z0
     *       x1 y1 z1
     *       x2 y2 x2
     *       x3 y3 x3
     *       x4 y4 x4
     *       x5 y5 x5
     *       x6 y6 x6
     *        0  0  0
     *
     *   if your structure has a number of atoms that is not a multiple of 4, you must
     *   pad it out to a multiple of 4 using zeros (using anything other than zero will
     *   make the calculation go wrong).
     *
     * On the other hand, when this file is compiled without -DALIGNED, then 
     * there are no 16 byte alignment or dummy atom requirements, and the
     * "npaddedatoms" argument is ignored.
     *
     *   arguments:
     *       nrealatoms:   the *actual* number of atoms in the structure
     *
     *       npaddedatoms: the number of atoms in the structure including padding atoms;
     *                     should equal nrealatoms rounded up to the next multiple of 4.
     *                     THIS IS ONLY USED IF COMPILED WITH -DALIGNED
     *
     *       a:            pointer to start of first structure (A). should be aligned to
     *                     a 16-byte boundary if compiled with -DALIGNED
     *
     *       b:            pointer to start of second structure (B). should be aligned to
     *                     a 16-byte boundary if compiled with -DALIGNED
     *
     *       G_a:          trace of A'A
     *
     *       G_b:          trace of B'B
     *
     *       computeRot:   if 0, the rotation matrix will not be computed. Otherwise, it
     *                     will be computed and stored in rot
     *
     *       rot:          output variable where the rotation matrix will be stored,
     *                     if computeRot != 0.
     */
    int niters, k;
#ifndef ALIGNED
    static const int masks[4][4] = {
        {1, 1, 1, 1},
        {1, 0, 0, 0},
        {1, 1, 0, 0},
        {1, 1, 1, 0}
    };
    int const *mask;
#endif
    /* Will have 3 garbage elements at the end */
    _ALIGNED(16) float M[12];
    __m128 xx,xy,xz,yx,yy,yz,zx,zy,zz;
    __m128 ax,ay,az,bx,by,bz;
    __m128 t0,t1,t2;

    if (a==b && G_a==G_b) {
        if (computeRot) {
            rot[0] = rot[4] = rot[8] = 1.0;
            rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;
        }
        return 0.0f;
    }

#ifdef ALIGNED
    niters = npaddedatoms >> 2;
    /* npaddedatoms must be a multiple of 4 */
    assert(npaddedatoms % 4 == 0);
#else
    niters = (nrealatoms + 4-1) / 4;
    mask = masks[nrealatoms%4];
#endif
    
    xx = xy = xz = yx = yy = yz = zx = zy = zz = _mm_setzero_ps();
    for (k = 0; k < niters; k++) {
#ifdef ALIGNED
        aos_deinterleaved_load(b,&bx,&by,&bz);
        aos_deinterleaved_load(a,&ax,&ay,&az);
#else 
        if (k == niters - 1) {
            /* x  y  z  */
            /* 0  1  2  */
            /* 3  4  5  */
            /* 6  7  8  */
            /* 9  10 11 */
            ax = _mm_set_ps(mask[0] ? a[0] : 0, mask[1] ? a[3] : 0, mask[2] ? a[6] : 0, mask[3] ? a[9] : 0);
            ay = _mm_set_ps(mask[0] ? a[1] : 0, mask[1] ? a[4] : 0, mask[2] ? a[7] : 0, mask[3] ? a[10] : 0);
            az = _mm_set_ps(mask[0] ? a[2] : 0, mask[1] ? a[5] : 0, mask[2] ? a[8] : 0, mask[3] ? a[11] : 0);

            bx = _mm_set_ps(mask[0] ? b[0] : 0, mask[1] ? b[3] : 0, mask[2] ? b[6] : 0, mask[3] ? b[9] : 0);
            by = _mm_set_ps(mask[0] ? b[1] : 0, mask[1] ? b[4] : 0, mask[2] ? b[7] : 0, mask[3] ? b[10] : 0);
            bz = _mm_set_ps(mask[0] ? b[2] : 0, mask[1] ? b[5] : 0, mask[2] ? b[8] : 0, mask[3] ? b[11] : 0);
        }
        else {
            aos_deinterleaved_loadu(b,&bx,&by,&bz);
            aos_deinterleaved_loadu(a,&ax,&ay,&az);
        }
#endif

        t0 = bx;
        t1 = by;
        t2 = bz;
        t0 = _mm_mul_ps(t0,ax);
        t1 = _mm_mul_ps(t1,ax);
        t2 = _mm_mul_ps(t2,ax);
        xx = _mm_add_ps(xx,t0);
        xy = _mm_add_ps(xy,t1);
        xz = _mm_add_ps(xz,t2);

        t0 = bx;
        t1 = by;
        t2 = bz;
        t0 = _mm_mul_ps(t0,ay);
        t1 = _mm_mul_ps(t1,ay);
        t2 = _mm_mul_ps(t2,ay);
        yx = _mm_add_ps(yx,t0);
        yy = _mm_add_ps(yy,t1);
        yz = _mm_add_ps(yz,t2);

        bx = _mm_mul_ps(bx,az);
        by = _mm_mul_ps(by,az);
        bz = _mm_mul_ps(bz,az);
        zx = _mm_add_ps(zx,bx);
        zy = _mm_add_ps(zy,by);
        zz = _mm_add_ps(zz,bz);

        a += 12;
        b += 12;
    }
    REDUCTION_EPILOGUE(xx, xy, xz, yx, yy, yz, zx, zy, zz, t0, t1, t2);

    _mm_store_ps(M  , xx);
    _mm_store_ps(M+4, yy);
    _mm_store_ps(M+8, zz);
    return msdFromMandG(M, G_a, G_b, nrealatoms, computeRot, rot);
}

