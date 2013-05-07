// This file is part of MSMBuilder.
//
// Copyright 2011 Stanford University
//
// MSMBuilder is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//
//=============================================================================================
// Calculation of RMSD by a the quaternion-based characteristic polynomial (QCP) algorithm of Theobald [1].
// 
// [1] Theobald DL. Rapid calculation of RMSDs using a quaternion-based characteristic polynomial. 
//     Acta Cryst., A61:478, 2005.  doi:10.1107/50108767305015266
//
// Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
// Contributions in 2010 from:
//      Kyle Beauchamp(kyleb@stanford.edu)
//      Peter Kasson (kasson@stanford.edu)
//      Kai Kohlhoff (kohlhoff@stanford.edu)
//      Imran Haque  (ihaque@cs.stanford.edu)
//=============================================================================================

#include <math.h>
#include "theobald_rmsd.h"
#include <stdio.h>
#include <xmmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif

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
				 int *nr12, int *nr34,double d0,double d1,double d2, double d3, double d4)
{
  double a3 = d3/d4;
  double a2 = d2/d4;
  double a1 = d1/d4;
  double a0 = d0/d4;

  double au2 = -a2;
  double au1 = (a1*a3 - 4.0*a0) ;
  double au0 = 4.0*a0*a2 - a1*a1 - a0*a3*a3;

  double x1, x2, x3;
  int nr = solve_cubic_equation(1.0, au2, au1, au0, &x1, &x2, &x3);

  double u1;
  if (nr==1) u1 = x1;
  else u1 = (x1>x3) ? x1 : x3;

  double R2 = 0.25*a3*a3 + u1 - a2;
  double R = (R2>0.0) ? sqrt(R2) : 0.0;

  double D2, E2;
  if (R != 0.0)
    {
      double foo1 = 0.75*a3*a3 - R2 - 2.0*a2;
      double foo2 = 0.25*(4.0*a3*a2 - 8.0*a1 - a3*a3*a3) / R;
      D2 = foo1 + foo2;
      E2 = foo1 - foo2;
    }
  else
    {
      double foo1 = 0.75*a3*a3 - 2.0*a2;
      double foo2 = 2.0 * sqrt(u1*u1 - 4.0*a0);
      D2 = foo1 + foo2;
      E2 = foo1 - foo2;
    }

  if (D2 >= 0.0)
    {
      double D = sqrt(D2);
      *r1 = -0.25*a3 + 0.5*R - 0.5*D;
      *r2 = -0.25*a3 + 0.5*R + 0.5*D;
      *nr12 = 2;
    }
  else
    {
      *r1 = *r2 = -0.25*a3 + 0.5*R;
      *nr12 = 0;
    }

  if (E2 >= 0.0)
    {
      double E = sqrt(E2);
      *r3 = -0.25*a3 - 0.5*R - 0.5*E;
      *r4 = -0.25*a3 - 0.5*R + 0.5*E;
      *nr34 = 2;
    }
  else
    {
      *r3 = *r4 = -0.25*a3 - 0.5*R;
      *nr34 = 0;
    }
  return *nr12 + *nr34;
}

int solve_cubic_equation(double  c3, double  c2,  double c1, double c0,
                         double *x1, double *x2, double *x3)
{
  double a2 = c2/c3;
  double a1 = c1/c3;
  double a0 = c0/c3;

  double q = a1/3.0 - a2*a2/9.0;
  double r = (a1*a2 - 3.0*a0)/6.0 - a2*a2*a2 / 27.0;
  double delta = q*q*q + r*r;

  if (delta>0.0)
    {
      double s1 = r + sqrt(delta);
      s1 = (s1>=0.0) ? pow(s1,1./3.) : -pow(-s1,1./3.);

      double s2 = r - sqrt(delta);
      s2 = (s2>=0.0) ? pow(s2,1./3.) : -pow(-s2,1./3.);

      *x1 = (s1+s2) - a2/3.0;
      *x2 = *x3 = -0.5 * (s1+s2) - a2/3.0;

      return 1;
    }
  else if (delta < 0.0)
    {
      double theta = acos(r/sqrt(-q*q*q)) / 3.0;
      double costh = cos(theta);
      double sinth = sin(theta);
      double sq = sqrt(-q);

      *x1 = 2.0*sq*costh - a2/3.0;
      *x2 = -sq*costh - a2/3.0 - sqrt(3.) * sq * sinth;
      *x3 = -sq*costh - a2/3.0 + sqrt(3.) * sq * sinth;

      return 3;
    }
  else
    {
      double s = (r>=0.0) ? pow(r,1./3.) : -pow(-r,1./3.);
      *x1 = 2.0*s - a2/3.0;
      *x2 = *x3 = -s - a2/3.0;

      return 3;
    }
}



float DirectSolve(float lambda, float C_0, float C_1, float C_2)
{
  double result;
  double r1,r2,r3,r4;
  int nr1,nr2;
  quartic_equation_solve_exact(&r1,&r2,&r3,&r4,&nr1,&nr2,(double )C_0,(double)C_1,(double)C_2,0.0,1.0);
  result=max(r1,r2);
  result=max(result,r3);
  result=max(result,r4);
  
  return(result);
}

float NewtonSolve(float lambda, float C_0, float C_1, float C_2)
{
  unsigned int maxits = 500;
  float tolerance = 1.0e-6f;
  float lambda_old,lambda2;
  float a,b;

  for (int i = 0; i < maxits; i++)
    {     
        lambda_old = lambda;
        lambda2 = lambda_old * lambda_old;
        b = (lambda2 + C_2) * lambda_old;
        a = b + C_1;
        lambda = lambda_old - (a * lambda_old + C_0) / (2.0f * lambda2 * lambda_old + b + a);
        if (fabsf(lambda - lambda_old) < fabsf(tolerance * lambda)) break;
    }
  if (fabsf(lambda - lambda_old) >= fabsf(100*tolerance * lambda))    
    {
      printf("RMSD Warning: No convergence after %d iterations: Lambda,Lambda0,Diff,Allowed = %f, %f, %f, %f \n",maxits,lambda, lambda_old, fabsf(lambda - lambda_old), fabsf(tolerance * lambda) );
    }
    
  return(lambda);
}

float rmsd2FromMandG(const float M[9],const float G_x,const float G_y,const int numAtoms) 
{
        // DEV: still too many glob mem accesses, can reduce this
    const int m = 3;
    float k00 =  M[0+0*m ] + M[1+1*m] + M[2+2*m];       // [0, 0]
    float k01 =  M[1+2*m ] - M[2+1*m];                  // [0, 1]
    float k02 =  M[2+0*m ] - M[0+2*m];                  // [0, 2]
    float k03 =  M[0+1*m ] - M[1+0*m];                  // [0, 3]
    float k11 =  M[0+0*m ] - M[1+1*m] - M[2+2*m];       // [1, 1]
    float k12 =  M[0+1*m ] + M[1+0*m];                  // [1, 2]
    float k13 =  M[2+0*m ] + M[0+2*m];                  // [1, 3]
    float k22 = -M[0+0*m ] + M[1+1*m] - M[2+2*m];       // [2, 2]
    float k23 =  M[1+2*m ] + M[2+1*m];                  // [2, 3]
    float k33 = -M[0+0*m ] - M[1+1*m] + M[2+2*m];       // [3, 3]


    // float C_4 = 1.0, C_3 = 0.0;
    float detM = 0.0f, detK = 0.0f;
    float C_2, C_1, C_0;

    float lambda = (G_x + G_y) / 2.0f;

    unsigned int i;
    float rmsd2,ls_rmsd2;
    C_2 = 0.0f;
    for (int i = 0; i < m * m; i++)
    {
        C_2 += M[i] * M[i];
    }
    C_2 *= -2.0f;

    // get determinante M
    // could use rule of Sarrus, but better:
    // computationally more efficient with Laplace expansion
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

    lambda=DirectSolve(lambda, C_0,C_1,C_2);

    rmsd2 = (G_x + G_y - 2.0f * lambda) / numAtoms;
    ls_rmsd2 = 0.0f;
    if (rmsd2 > 0.0f) ls_rmsd2 = rmsd2;

    return ls_rmsd2;
    
}

real ls_rmsd2_aligned_T_g(const int nrealatoms, const int npaddedatoms, const int rowstride,
                          const float* aT, const float* bT, const real G_a, const real G_b)
{
    /* Structure setup for this function:
     *
     *   structures are stored axis major, possibly with extra padding to ensure you
     *   meet two constraints:
     *       - the number of elements in a row must be a multiple of 4
     *       - the first element in each row must be aligned to a 16 byte boundary
     *
     *   note that if you meet the second condition for the first row, and meet the
     *   first condition, the alignment will automatically be satisfied for every row.
     *
     *   the layout in memory for a structure of 7 atoms would look like this:
     *
     *       x0 x1 x2 x3 x4 x5 x6 0
     *       y0 y1 y2 y3 y4 y5 y6 0
     *       z0 z1 z2 z3 z4 z5 z6 0
     *
     *   if your structure has a number of atoms that is not a multiple of 4, you must
     *   pad it out to a multiple of 4 using zeros (using anything other than zero will
     *   make the calculation go wrong).
     *
     *   arguments:
     *       nrealatoms:   the *actual* number of atoms in the structure
     *
     *       npaddedatoms: the number of atoms in the structure including padding atoms;
     *                     should equal nrealatoms rounded up to the next multiple of 4
     *
     *       rowstride:    the offset in elements between rows in the arrays. will prob
     *                     be equal to npaddedatoms, but you might use something else if
     *                     (for example) you were subsetting the structure
     *
     *       aT:           pointer to start of first structure (A). should be aligned to
     *                     a 16-byte boundary
     *
     *       bT:           pointer to start of second structure (B). should be aligned to
     *                     a 16-byte boundary
     *
     *       G_a:          trace of A'A
     *
     *       G_b:          trace of B'B
     */

            
            
	int nIndex;
    // Will have 3 garbage elements at the end
    float M[12] __attribute__ ((aligned (16)));

    const float* aTx = aT;
    const float* aTy = aT+rowstride;
    const float* aTz = aT+2*rowstride;
    const float* bTx = bT;
    const float* bTy = bT+rowstride;
    const float* bTz = bT+2*rowstride;


    // npaddedatoms must be a multiple of 4
    int niters = npaddedatoms >> 2;
    __m128 xx,xy,xz,yx,yy,yz,zx,zy,zz;
    __m128 ax,ay,az,b;
    __m128 t0,t1,t2;
    // Prologue
    xx = _mm_xor_ps(xx,xx);
    xy = _mm_xor_ps(xy,xy);
    xz = _mm_xor_ps(xz,xz);
    yx = _mm_xor_ps(yx,yx);
    yy = _mm_xor_ps(yy,yy);
    yz = _mm_xor_ps(yz,yz);
    zx = _mm_xor_ps(zx,zx);
    zy = _mm_xor_ps(zy,zy);
    zz = _mm_xor_ps(zz,zz);
    for (int k = 0; k < niters; k++) {
        ax = _mm_load_ps(aTx);
        ay = _mm_load_ps(aTy);
        az = _mm_load_ps(aTz);

        b = _mm_load_ps(bTx);
        t0 = ax;
        t1 = ay;
        t2 = az;

        t0 = _mm_mul_ps(t0,b);
        t1 = _mm_mul_ps(t1,b);
        t2 = _mm_mul_ps(t2,b);

        xx = _mm_add_ps(xx,t0);
        yx = _mm_add_ps(yx,t1);
        zx = _mm_add_ps(zx,t2);

        b = _mm_load_ps(bTy);
        t0 = ax;
        t1 = ay;
        t2 = az;

        t0 = _mm_mul_ps(t0,b);
        t1 = _mm_mul_ps(t1,b);
        t2 = _mm_mul_ps(t2,b);

        xy = _mm_add_ps(xy,t0);
        yy = _mm_add_ps(yy,t1);
        zy = _mm_add_ps(zy,t2);

        b = _mm_load_ps(bTz);

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
    // Epilogue - reduce 4 wide vectors to one wide
   /*xmm07 = xx0 xx1 xx2 xx3
     xmm08 = xy0 xy1 xy2 xy3
     xmm09 = xz0 xz1 xz2 xz3
     xmm10 = yx0 yx1 yx2 yx3
     xmm11 = yy0 yy1 yy2 yy3
     xmm12 = yz0 yz1 yz2 yz3
     xmm13 = zx0 zx1 zx2 zx3
     xmm14 = zy0 zy1 zy2 zy3
     xmm15 = zz0 zz1 zz2 zz3
     
     haddps xmm07 xmm08
         xmm07 = xx0+1 xx2+3 xy0+1 xy2+3
     haddps xmm09 xmm10
         xmm09 = xz0+1 xz2+3 yx0+1 yx2+3
     haddps xmm11 xmm12
         xmm11 = yy0+1 yy2+3 yz0+1 yz2+3
     haddps xmm13 xmm14
         xmm13 = zx0+1 zx2+3 zy0+1 zy2+3
     haddps xmm15 xmm14
         xmm15 = zz0+1 zz2+3 zy0+1 zy2+3
     
     haddps xmm07 xmm09
         xmm07 = xx0123 xy0123 xz0123 yx0123
     haddps xmm11 xmm13
         xmm11 = yy0123 yz0123 zx0123 zy0123
     haddps xmm15 xmm09
         xmm15 = zz0123 zy0123 xz0123 yx0123*/ 
    #ifdef __SSE3__
    xx = _mm_hadd_ps(xx,xy);
    xz = _mm_hadd_ps(xz,yx);
    yy = _mm_hadd_ps(yy,yz);
    zx = _mm_hadd_ps(zx,zy);
    zz = _mm_hadd_ps(zz,zy);

    xx = _mm_hadd_ps(xx,xz);
    yy = _mm_hadd_ps(yy,zx);
    zz = _mm_hadd_ps(zz,xz);
    #else
    // Emulate horizontal adds using UNPCKLPS/UNPCKHPS
    t0 = xx;
    t1 = xx;
    t0 = _mm_unpacklo_ps(t0,xz);
        // = xx0 xz0 xx1 xz1
    t1 = _mm_unpackhi_ps(t1,xz);
        // = xx2 xz2 xx3 xz3
    t0 = _mm_add_ps(t0,t1);
        // = xx02 xz02 xx13 xz13

    t1 = xy;
    t2 = xy;
    t1 = _mm_unpacklo_ps(t1,yx);
        // = xy0 yx0 xy1 yx1
    t2 = _mm_unpackhi_ps(t2,yx);
        // = xy2 yx2 xy3 yx3
    t1 = _mm_add_ps(t1,t2);
        // = xy02 yx02 xy13 yx13

    xx = t0;
    xx = _mm_unpacklo_ps(xx,t1);
        // = xx02 xy02 xz02 yx02
    t0 = _mm_unpackhi_ps(t0,t1);
        // = xx13 xy13 xz13 yx13
    xx = _mm_add_ps(xx,t0);
        // = xx0123 xy0123 xz0123 yx0123

    t0 = yy;
    t1 = yy;
    t0 = _mm_unpacklo_ps(t0,zx);
        // = yy0 zx0 yy1 zx1
    t1 = _mm_unpackhi_ps(t1,zx);
        // = yy2 zx2 yy3 zx3
    t0 = _mm_add_ps(t0,t1);
        // = yy02 zx02 yy13 zx13

    t1 = yz;
    t2 = yz;
    t1 = _mm_unpacklo_ps(t1,zy);
        // = yz0 zy0 yz1 zy1
    t2 = _mm_unpackhi_ps(t2,zy);
        // = yz2 zy2 yz3 zy3
    t1 = _mm_add_ps(t1,t2);
        // = yz02 zy02 yz13 zy13

    yy = t0;
    yy = _mm_unpacklo_ps(yy,t1);
        // = yy02 yz02 zx02 zy02
    t0 = _mm_unpackhi_ps(t0,t1);
        // = yy13 yz13 zx13 zy13
    yy = _mm_add_ps(yy,t0);
        // = yy0123 yz0123 zx0123 zy0123

    t1 = _mm_movehl_ps(t1,zz);
        // = zz2 zz3 - -
    zz = _mm_add_ps(zz,t1);
        // = zz02 zz13 - -
    t1 = _mm_shuffle_ps(zz,zz,_MM_SHUFFLE(1,1,1,1));
        // = zz13 zz13 zz13 zz13
    zz = _mm_add_ps(zz,t1);
        // = zz0123 zz1133 - -
    #endif

    _mm_store_ps(M  , xx);
    _mm_store_ps(M+4, yy);
    _mm_store_ps(M+8, zz);

    return rmsd2FromMandG(M,G_a,G_b,nrealatoms);
}

