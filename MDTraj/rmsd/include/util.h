#ifndef __UTIL_H__
#define __UTIL_H__
/*
    Utility functions for IRMSD matrix multiply routines.

    Common intrinsic subroutines shared among files.
*/
#include "msvccompat.h"
#include <xmmintrin.h>
#ifdef __SSSE3__
#include <pmmintrin.h>
#endif
#include "sse_swizzle.h"

static INLINE void aos_interleaved_storeu(float* p, __m128 x, __m128 y, __m128 z)
{
    /* Store the coordinates of 4 atoms currently in struct-of-arrays format
       in an array-of-struct memory layout, using unaligned memory

    p points to memory that is supposed to be layed out in atom major order, like

    x1 y1 z1
    x2 y2 z2
    x3 y3 z3
    x4 y4 z4
    ...

    And x, y, and z are SSE float4s that contain the x, y and z coordinates of
    four atoms. In order to actually store the 12 floats in memory, you need to
    do some shuffles first.

    This code is directly translated from the Intel link below:
    http://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx
    */
    __m128 x0x2y0y2 = _mm_shuffle_ps(x, y, _MM_SHUFFLE(2,0,2,0));
    __m128 y1y3z1z3 = _mm_shuffle_ps(y, z, _MM_SHUFFLE(3,1,3,1));
    __m128 z0z2x1x3 = _mm_shuffle_ps(z, x, _MM_SHUFFLE(3,1,2,0));

    __m128 rx0y0z0x1 = _mm_shuffle_ps(x0x2y0y2, z0z2x1x3, _MM_SHUFFLE(2,0,2,0));
    __m128 ry1z2x2y2 = _mm_shuffle_ps(y1y3z1z3, x0x2y0y2, _MM_SHUFFLE(3,1,2,0));
    __m128 rz2x3y3z3 = _mm_shuffle_ps(z0z2x1x3, y1y3z1z3, _MM_SHUFFLE(3,1,3,1));

    _mm_storeu_ps(p+0, rx0y0z0x1);
    _mm_storeu_ps(p+4, ry1z2x2y2);
    _mm_storeu_ps(p+8, rz2x3y3z3);
}

static INLINE void aos_interleaved_store(float* p, __m128 x, __m128 y, __m128 z)
{
     /* Store the coordinates of 4 atoms currently in struct-of-arrays format
        in an array-of-struct memory layout, using aligned memory

     p points to memory that is supposed to be layed out in atom major order, like

     x1 y1 z1
     x2 y2 z2
     x3 y3 z3
     x4 y4 z4
     ...

     And x, y, and z are SSE float4s that contain the x, y and z coordinates of
     four atoms. In order to actually store the 12 floats in memory, you need to
     do some shuffles first.

     This code is directly translated from the Intel link below:
     http://software.intel.com/en-us/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx
     */
    __m128 x0x2y0y2 = _mm_shuffle_ps(x, y, _MM_SHUFFLE(2,0,2,0));
    __m128 y1y3z1z3 = _mm_shuffle_ps(y, z, _MM_SHUFFLE(3,1,3,1));
    __m128 z0z2x1x3 = _mm_shuffle_ps(z, x, _MM_SHUFFLE(3,1,2,0));

    __m128 rx0y0z0x1 = _mm_shuffle_ps(x0x2y0y2, z0z2x1x3, _MM_SHUFFLE(2,0,2,0));
    __m128 ry1z2x2y2 = _mm_shuffle_ps(y1y3z1z3, x0x2y0y2, _MM_SHUFFLE(3,1,2,0));
    __m128 rz2x3y3z3 = _mm_shuffle_ps(z0z2x1x3, y1y3z1z3, _MM_SHUFFLE(3,1,3,1));

    _mm_store_ps(p+0, rx0y0z0x1);
    _mm_store_ps(p+4, ry1z2x2y2);
    _mm_store_ps(p+8, rz2x3y3z3);
}


static INLINE void aos_deinterleaved_load(const float* S, __m128* x, __m128* y, __m128* z)
{
    /* Deinterleaved load of x, y, z coordinates of 4 atoms in atom-major
       ordering using aligned loads.
    */

  __m128 t1, t2;
  *x  = _mm_load_ps(S);
  *y  = _mm_load_ps(S+4);
  t1 = _mm_load_ps(S+8);
  *z  = *x;
  t2 = *y;

  t2 = _mm_shuffle_ps_yzyz(t2,t1);
  *z  = _mm_shuffle_ps_yzxw(*z,t2);
  *x  = _mm_shuffle_ps_xwyz(*x,t2);
  *y  = _mm_shuffle_ps_xwxw(*y,*z);
  *z  = _mm_shuffle_ps_yzxw(*z,t1);
  *y  = _mm_swizzle_ps_zxyw(*y);
  return;
}

static INLINE void aos_deinterleaved_loadu(const float* S, __m128* x, __m128* y, __m128* z)
{
    /* Deinterleaved load of x, y, z coordinates of 4 atoms in atom-major
       ordering using unaligned loads
    */

  __m128 t1, t2;
  *x  = _mm_loadu_ps(S);
  *y  = _mm_loadu_ps(S+4);
  t1 = _mm_loadu_ps(S+8);
  *z  = *x;
  t2 = *y;

  t2 = _mm_shuffle_ps_yzyz(t2,t1);
  *z  = _mm_shuffle_ps_yzxw(*z,t2);
  *x  = _mm_shuffle_ps_xwyz(*x,t2);
  *y  = _mm_shuffle_ps_xwxw(*y,*z);
  *z  = _mm_shuffle_ps_yzxw(*z,t1);
  *y  = _mm_swizzle_ps_zxyw(*y);
  return;
}

/* Epilogue - reduce 4 wide vectors to one wide */
#ifdef __SSE3__
/* Use SSE3 horizontal add to do the reduction */
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
      xmm15 = zz0123 zy0123 xz0123 yx0123
*/
#define REDUCTION_EPILOGUE(xx, xy, xz, yx, yy, yz, zx, zy, zz, t0, t1, t2) \
  xx = _mm_hadd_ps(xx,xy); /* xx = xx01 xx23 xy01 xy23 */\
  xz = _mm_hadd_ps(xz,yx); /* xz = xz01 xz23 yx01 yx23 */\
  yy = _mm_hadd_ps(yy,yz); /* yy = yy01 yy23 yz01 yz23 */\
  zx = _mm_hadd_ps(zx,zy); /* zx = zx01 zx23 zy01 zy23 */\
  zz = _mm_hadd_ps(zz,zy); /* zz = zz01 zz23 zy01 zy23 */\
  xx = _mm_hadd_ps(xx,xz); /* xx = xx0123 xy0123 xz0123 yx0123 */\
  yy = _mm_hadd_ps(yy,zx); /* yy = yy0123 yz0123 zx0123 zy0123 */\
  zz = _mm_hadd_ps(zz,xz); /* zz = zz0123 zy0123 xz0123 yx0123 */\

#else
/* Emulate horizontal adds using SSE2 UNPCKLPS/UNPCKHPS */
#define REDUCTION_EPILOGUE(xx, xy, xz, yx, yy, yz, zx, zy, zz, t0, t1, t2) \
  t0 = xx; \
  t1 = xx; \
  t0 = _mm_unpacklo_ps(t0,xz); \
  /* = xx0 xz0 xx1 xz1 */ \
  t1 = _mm_unpackhi_ps(t1,xz); \
  /* = xx2 xz2 xx3 xz3 */ \
  t0 = _mm_add_ps(t0,t1); \
  /* = xx02 xz02 xx13 xz13 */ \
    \
  t1 = xy; \
  t2 = xy; \
  t1 = _mm_unpacklo_ps(t1,yx); \
  /* = xy0 yx0 xy1 yx1 */ \
  t2 = _mm_unpackhi_ps(t2,yx); \
  /* = xy2 yx2 xy3 yx3 */ \
  t1 = _mm_add_ps(t1,t2); \
  /* = xy02 yx02 xy13 yx13 */ \
    \
  xx = t0; \
  xx = _mm_unpacklo_ps(xx,t1); \
  /* = xx02 xy02 xz02 yx02 */ \
  t0 = _mm_unpackhi_ps(t0,t1); \
  /* = xx13 xy13 xz13 yx13 */ \
  xx = _mm_add_ps(xx,t0); \
  /* = xx0123 xy0123 xz0123 yx0123 */ \
    \
  t0 = yy; \
  t1 = yy; \
  t0 = _mm_unpacklo_ps(t0,zx); \
  /* = yy0 zx0 yy1 zx1 */ \
  t1 = _mm_unpackhi_ps(t1,zx); \
  /* = yy2 zx2 yy3 zx3 */ \
  t0 = _mm_add_ps(t0,t1); \
  /* = yy02 zx02 yy13 zx13 */ \
    \
  t1 = yz; \
  t2 = yz; \
  t1 = _mm_unpacklo_ps(t1,zy); \
  /* = yz0 zy0 yz1 zy1 */ \
  t2 = _mm_unpackhi_ps(t2,zy); \
  /* = yz2 zy2 yz3 zy3 */ \
  t1 = _mm_add_ps(t1,t2); \
  /* = yz02 zy02 yz13 zy13 */ \
    \
  yy = t0; \
  yy = _mm_unpacklo_ps(yy,t1); \
  /* = yy02 yz02 zx02 zy02 */ \
  t0 = _mm_unpackhi_ps(t0,t1); \
  /* = yy13 yz13 zx13 zy13 */ \
  yy = _mm_add_ps(yy,t0); \
  /* = yy0123 yz0123 zx0123 zy0123 */ \
    \
  t1 = _mm_movehl_ps(t1,zz); \
  /* = zz2 zz3 - - */ \
  zz = _mm_add_ps(zz,t1); \
  /* = zz02 zz13 - - */ \
  t1 = _mm_shuffle_ps(zz,zz,_MM_SHUFFLE(1,1,1,1)); \
  /* = zz13 zz13 zz13 zz13 */ \
  zz = _mm_add_ps(zz,t1); \
  /* = zz0123 zz1133 - - */
#endif

#endif
