/*/////////////////////////////////////////////////////////////////////////////
// MDTraj: A Python Library for Loading, Saving, and Manipulating
//         Molecular Dynamics Trajectories.
// Copyright 2012-2013 Stanford University and the Authors
//
// Authors: Imran S. Haque
// Contributors: John D. Chodera, Kyle Beauchamp, Peter Kasson,
//               Kai Kohlhoff, Jonathan Zrake, Robert McGibbon
//
// MDTraj is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Calculation of RMSD by a the quaternion-based characteristic polynomial
// (QCP) algorithm of Theobald [1].
//
// [1] Theobald DL. Rapid calculation of RMSDs using a quaternion-based
//     characteristic polynomial. Acta Cryst., A61:478, 2005.
//     doi:10.1107/50108767305015266
///////////////////////////////////////////////////////////////////////////////*/

#include "msvccompat.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <xmmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#include "theobald_rmsd.h"
#include "util.h"
#include "stdio.h"

/* 
Using #define ALIGNED enables aligned loads. Slight speed boost, but
requires that input data be 16-byte aligned.
*/
/* #define ALIGNED */

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif


/*------------------------------------------------------------------------------
 * The quartic and cubic functions are taken from:
 * FILE: quartic.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *         Adapted from the nvwa code by Weiqun Zhang
 * Modified by KAB 2011
 * GPLv2 / LGPL exemption from Jonathan Zrake, Aug. 2, 2011
 * Original code from http://code.google.com/p/python-mhd/
 *------------------------------------------------------------------------------
 */

int quartic_equation_solve_exact(double *r1, double *r2, double *r3, double *r4,
                                 int *nr12, int *nr34,double d0,double d1,double d2, double d3, double d4) {
    int nr;
    double a0, a1, a2, a3, au0, au1, au2;
    double x1, x2, x3, u1, R, R2, D, D2, E, E2, foo1, foo2;

    a3 = d3/d4;
    a3 = d3/d4;
    a2 = d2/d4;
    a1 = d1/d4;
    a0 = d0/d4;

    au2 = -a2;
    au1 = (a1*a3 - 4.0*a0) ;
    au0 = 4.0*a0*a2 - a1*a1 - a0*a3*a3;

    nr = solve_cubic_equation(1.0, au2, au1, au0, &x1, &x2, &x3);

    if (nr==1) u1 = x1;
    else u1 = (x1>x3) ? x1 : x3;

    R2 = 0.25*a3*a3 + u1 - a2;
    R = (R2>0.0) ? sqrt(R2) : 0.0;

    if (R != 0.0) {
        foo1 = 0.75*a3*a3 - R2 - 2.0*a2;
        foo2 = 0.25*(4.0*a3*a2 - 8.0*a1 - a3*a3*a3) / R;
        D2 = foo1 + foo2;
        E2 = foo1 - foo2;
    } else {
        foo1 = 0.75*a3*a3 - 2.0*a2;
        foo2 = 2.0 * sqrt(u1*u1 - 4.0*a0);
        D2 = foo1 + foo2;
        E2 = foo1 - foo2;
    }

    if (D2 >= 0.0) {
        D = sqrt(D2);
        *r1 = -0.25*a3 + 0.5*R - 0.5*D;
        *r2 = -0.25*a3 + 0.5*R + 0.5*D;
        *nr12 = 2;
    } else {
        *r1 = *r2 = -0.25*a3 + 0.5*R;
        *nr12 = 0;
    }

    if (E2 >= 0.0) {
        E = sqrt(E2);
        *r3 = -0.25*a3 - 0.5*R - 0.5*E;
        *r4 = -0.25*a3 - 0.5*R + 0.5*E;
        *nr34 = 2;
    } else {
        *r3 = *r4 = -0.25*a3 - 0.5*R;
        *nr34 = 0;
    }
    return *nr12 + *nr34;
}

int solve_cubic_equation(double  c3, double  c2,  double c1, double c0,
                         double *x1, double *x2, double *x3) {
    double s1, s2;
    double a2 = c2/c3;
    double a1 = c1/c3;
    double a0 = c0/c3;

    double q = a1/3.0 - a2*a2/9.0;
    double r = (a1*a2 - 3.0*a0)/6.0 - a2*a2*a2 / 27.0;
    double delta = q*q*q + r*r;

    if (delta>0.0) {
        s1 = r + sqrt(delta);
        s1 = (s1>=0.0) ? pow(s1,1./3.) : -pow(-s1,1./3.);

        s2 = r - sqrt(delta);
        s2 = (s2>=0.0) ? pow(s2,1./3.) : -pow(-s2,1./3.);

        *x1 = (s1+s2) - a2/3.0;
        *x2 = *x3 = -0.5 * (s1+s2) - a2/3.0;

        return 1;
    } else if (delta < 0.0) {
        double theta = acos(r/sqrt(-q*q*q)) / 3.0;
        double costh = cos(theta);
        double sinth = sin(theta);
        double sq = sqrt(-q);

        *x1 = 2.0*sq*costh - a2/3.0;
        *x2 = -sq*costh - a2/3.0 - sqrt(3.) * sq * sinth;
        *x3 = -sq*costh - a2/3.0 + sqrt(3.) * sq * sinth;

        return 3;
    } else {
        double s = (r>=0.0) ? pow(r,1./3.) : -pow(-r,1./3.);
        *x1 = 2.0*s - a2/3.0;
        *x2 = *x3 = -s - a2/3.0;

        return 3;
    }
}



float DirectSolve(float lambda, float C_0, float C_1, float C_2) {
    double result;
    double r1,r2,r3,r4;
    int nr1,nr2;
    quartic_equation_solve_exact(&r1,&r2,&r3,&r4,&nr1,&nr2,(double )C_0,(double)C_1,(double)C_2,0.0,1.0);
    result=max(r1,r2);
    result=max(result,r3);
    result=max(result,r4);

    return (float) result;
}

float NewtonSolve(float lambda, float C_0, float C_1, float C_2) {
    unsigned int i;
    unsigned int maxits = 500;
    float tolerance = 1.0e-6f;
    float lambda_old,lambda2;
    float a,b;

    for (i = 0; i < maxits; i++) {
        lambda_old = lambda;
        lambda2 = lambda_old * lambda_old;
        b = (lambda2 + C_2) * lambda_old;
        a = b + C_1;
        lambda = lambda_old - (a * lambda_old + C_0) / (2.0f * lambda2 * lambda_old + b + a);
        if (fabsf(lambda - lambda_old) < fabsf(tolerance * lambda)) break;
    }
    if (fabsf(lambda - lambda_old) >= fabsf(100*tolerance * lambda)) {
        printf("RMSD Warning: No convergence after %d iterations: Lambda,Lambda0,Diff,Allowed = %f, %f, %f, %f \n",maxits,lambda, lambda_old, fabsf(lambda - lambda_old), fabsf(tolerance * lambda) );
    }

    return(lambda);
}

float msdFromMandG(const float M[9], const float G_x, const float G_y,
                   const int numAtoms, int computeRot, float rot[9]) {
    /* Compute the RMSD and (optionally) the rotation matrix from M and G. Core routine of the
       theobald QCP method.

    when computeRot == 0, the rotation matrix will not be computed.
    otherwise, it will be computed and stored in `rot`.
    */
    int i;
    const int m = 3;
    float k00 =  M[0+0*m ] + M[1+1*m] + M[2+2*m];       /* [0, 0] */
    float k01 =  M[1+2*m ] - M[2+1*m];                  /* [0, 1] */
    float k02 =  M[2+0*m ] - M[0+2*m];                  /* [0, 2] */
    float k03 =  M[0+1*m ] - M[1+0*m];                  /* [0, 3] */
    float k11 =  M[0+0*m ] - M[1+1*m] - M[2+2*m];       /* [1, 1] */
    float k12 =  M[0+1*m ] + M[1+0*m];                  /* [1, 2] */
    float k13 =  M[2+0*m ] + M[0+2*m];                  /* [1, 3] */
    float k22 = -M[0+0*m ] + M[1+1*m] - M[2+2*m];       /* [2, 2] */
    float k23 =  M[1+2*m ] + M[2+1*m];                  /* [2, 3] */
    float k33 = -M[0+0*m ] - M[1+1*m] + M[2+2*m];       /* [3, 3] */

    /* declarations required for rotation code */
    float k2233_2323, k1233_1323, k1223_1322, k0223_0322, k0233_0323, k0213_0312;
    float q0, q1, q2, q3, qsqr, normq, a2, x2, y2, z2, xy, az, zx, ay, yz, ax;

    /* float C_4 = 1.0, C_3 = 0.0; */
    float detM = 0.0f, detK = 0.0f;
    float C_2, C_1, C_0;

    float lambda = (G_x + G_y) / 2.0f;

    float rmsd2,ls_rmsd2;
    C_2 = 0.0f;
    for (i = 0; i < m * m; i++) {
        C_2 += M[i] * M[i];
    }
    C_2 *= -2.0f;

    /* get det(M) */
    /* could use rule of Sarrus, but better: */
    /* computationally more efficient with Laplace expansion */
    detM = M[0] * (M[4] * M[8] - M[5] * M[7])
           + M[3] * (M[7] * M[2] - M[8] * M[1])
           + M[6] * (M[1] * M[5] - M[2] * M[4]);

    detK = k01*k01*k23*k23   - k22*k33*k01*k01   + 2*k33*k01*k02*k12
           - 2*k01*k02*k13*k23 - 2*k01*k03*k12*k23 + 2*k22*k01*k03*k13
           + k02*k02*k13*k13   - k11*k33*k02*k02   - 2*k02*k03*k12*k13
           + 2*k11*k02*k03*k23 + k03*k03*k12*k12   - k11*k22*k03*k03
           - k00*k33*k12*k12   + 2*k00*k12*k13*k23 - k00*k22*k13*k13
           - k00*k11*k23*k23   + k00*k11*k22*k33;


    C_1 = -8.0f * detM;
    C_0 = detK;

    lambda = DirectSolve(lambda, C_0,C_1,C_2);

    rmsd2 = (G_x + G_y - 2.0f * lambda) / numAtoms;
    ls_rmsd2 = 0.0f;
    if (rmsd2 > 0.0f) ls_rmsd2 = rmsd2;

    if (computeRot != 0) {
        /* compute the rotation matrix */
        k00 -= lambda;
        k11 -= lambda;
        k22 -= lambda;
        k33 -= lambda;

        k2233_2323 = k22*k33 - k23*k23;
        k1233_1323 = k12*k33 - k13*k23;
        k1223_1322 = k12*k23 - k13*k22;
        k0223_0322 = k02*k23 - k03*k22;
        k0233_0323 = k02*k33 - k03*k23;
        k0213_0312 = k02*k13 - k03*k12;

        q0 =  k11*k2233_2323 - k12*k1233_1323 + k13*k1223_1322;
        q1 = -k01*k2233_2323 + k12*k0233_0323 - k13*k0223_0322;
        q2 =  k01*k1233_1323 - k11*k0233_0323 + k13*k0213_0312;
        q3 = -k01*k1223_1322 + k11*k0223_0322 - k12*k0213_0312;
        qsqr = q0*q0 + q1*q1 + q2*q2 + q3*q3;

        if (qsqr < 1e-11f) {
            fprintf(stderr, "%s UNCONVERGED ROTATION MATRIX. RETURNING IDENTITY=%d\n", __FILE__, __LINE__);
            rot[0] = rot[4] = rot[8] = 1.0;
            rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;
        } else {
            normq = sqrt(qsqr);
            q0 /= normq;
            q1 /= normq;
            q2 /= normq;
            q3 /= normq;

            a2 = q0 * q0;
            x2 = q1 * q1;
            y2 = q2 * q2;
            z2 = q3 * q3;

            xy = q1 * q2;
            az = q0 * q3;
            zx = q3 * q1;
            ay = q0 * q2;
            yz = q2 * q3;
            ax = q0 * q1;

            rot[0] = a2 + x2 - y2 - z2;
            rot[3] = 2 * (xy + az);
            rot[6] = 2 * (zx - ay);
            rot[1] = 2 * (xy - az);
            rot[4] = a2 - x2 + y2 - z2;
            rot[7] = 2 * (yz + ax);
            rot[2] = 2 * (zx + ay);
            rot[5] = 2 * (yz - ax);
            rot[8] = a2 - x2 - y2 + z2;
        }
    }
    return ls_rmsd2;
}

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

