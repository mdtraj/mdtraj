#ifndef _SSE_TOOLS_H_
#define _SSE_TOOLS_H_

#include "msvccompat.h"
#include <stdio.h>
#include <pmmintrin.h>

/* Macros from http://sseplus.sourceforge.net/_s_s_e_plus__platform_8h-source.html */

#ifdef _MSC_VER
#define __CNST32I28I_( x ) \
    ((unsigned __int8)((x) & 0xFF)), ((unsigned __int8)(((x) >> 8) & 0xFF)), ((unsigned __int8)(((x) >> 16) & 0xFF)), ((unsigned __int8)(((x) >> 24) & 0xFF))

#define SSP_CONST_SETR_32I( a, b, c, d ) \
    { __CNST32I28I_((a)), __CNST32I28I_((b)), __CNST32I28I_((c)), __CNST32I28I_((d)) }

#define SSP_CONST_SET_32I( a, b, c, d ) \
    SSP_CONST_SETR_32I( (d), (c), (b), (a) )
#else
#define __CNST32TO64_( a, b ) \
        ( ((b)<<32) | ((a) & 0xFFFFFFFF) )

#define SSP_CONST_SETR_32I( a, b, c, d ) \
    { __CNST32TO64_( (unsigned long long)(a), (unsigned long long)(b) ), \
      __CNST32TO64_( (unsigned long long)(c), (unsigned long long)(d) ) }

#define SSP_CONST_SET_32I( a, b, c, d ) \
    SSP_CONST_SETR_32I( (d), (c), (b), (a) )
#endif


/****************************************************************************/
/* Utilities                                                                */
/****************************************************************************/

static INLINE __m128 _mm_hsum_ps(__m128 v) {
    v = _mm_hadd_ps(v, v);
    v = _mm_hadd_ps(v, v);
    return v;
}

static INLINE __m128 _mm_round_ps2(const __m128 a){
    /* http://dss.stephanierct.com/DevBlog/?p=8 */
    __m128 v0 = _mm_setzero_ps();             /* generate the highest value < 2 */
    __m128 v1 = _mm_cmpeq_ps(v0,v0);
    __m128i srli = _mm_srli_epi32( *(__m128i*)& v1, 2);
    __m128 vNearest2 = *(__m128*)& srli;
    __m128i i = _mm_cvttps_epi32(a);
    __m128 aTrunc = _mm_cvtepi32_ps(i);        /* truncate a */
    __m128 rmd = _mm_sub_ps(a, aTrunc);        /* get remainder */
    __m128 rmd2 = _mm_mul_ps( rmd, vNearest2); /* mul remainder by near 2 will yield the needed offset */
    __m128i rmd2i = _mm_cvttps_epi32(rmd2);    /* after being truncated of course */
    __m128 rmd2Trunc = _mm_cvtepi32_ps(rmd2i);
    __m128 r =_mm_add_ps(aTrunc, rmd2Trunc);
    return r;
}


static INLINE __m128 _mm_dp_ps2( __m128 a, __m128 b, const int mask ) {
    /*
     Copyright (c) 2006-2008 Advanced Micro Devices, Inc. All Rights Reserved.
     This software is subject to the Apache v2.0 License.
    */

    /* Shift mask multiply moves 0,1,2,3 bits to left, becomes MSB */
    const static __m128i mulShiftImm_0123 = SSP_CONST_SET_32I(0x010000, 0x020000, 0x040000, 0x080000);
    /* Shift mask multiply moves 4,5,6,7 bits to left, becomes MSB */
    const static __m128i mulShiftImm_4567 = SSP_CONST_SET_32I(0x100000, 0x200000, 0x400000, 0x800000);

    /* Begin mask preparation */
    __m128i mHi, mLo;
    mLo = _mm_set1_epi32(mask);    /* Load the mask into register */
    mLo = _mm_slli_si128(mLo, 3);  /* Shift into reach of the 16 bit multiply */
    mHi = _mm_mullo_epi16(mLo, mulShiftImm_0123);  /* Shift the bits */
    mLo = _mm_mullo_epi16(mLo, mulShiftImm_4567);  /* Shift the bits */
    mHi = _mm_cmplt_epi32(mHi, _mm_setzero_si128()); /* FFFFFFFF if bit set, 00000000 if not set */
    mLo = _mm_cmplt_epi32(mLo, _mm_setzero_si128()); /* FFFFFFFF if bit set, 00000000 if not set */
    /* End mask preparation - Mask bits 0-3 in mLo, 4-7 in mHi */
    a = _mm_and_ps(a, *(__m128*)& mHi);   /* Clear input using the high bits of the mask */
    a = _mm_mul_ps(a, b);
    a = _mm_hsum_ps(a);                  /* Horizontally add the 4 values */
    a = _mm_and_ps(a, *(__m128*)& mLo);  /* Clear output using low bits of the mask */
    return a;
}


static INLINE __m128 load_float3(const float* value) {
    /* Load (x,y,z) into a SSE register, leaving the last entry */
    /* set to zero. */
  __m128 x = _mm_load_ss(&value[0]);
  __m128 y = _mm_load_ss(&value[1]);
  __m128 z = _mm_load_ss(&value[2]);
  __m128 xy = _mm_movelh_ps(x, y);
  return _mm_shuffle_ps(xy, z, _MM_SHUFFLE(2, 0, 2, 0));
}

static INLINE int store_float3(float* loc, __m128 val) {
    /* Store the low three floats in an SSE register into */
    /* memory, at location loc[0], loc[1], loc[2]. The high */
    /* float is not touched. */
  _mm_store_ss(loc, val);
  _mm_store_ss(loc+1, _mm_shuffle_ps(val, val, _MM_SHUFFLE(1,1,1,1)));
  _mm_store_ss(loc+2, _mm_shuffle_ps(val, val, _MM_SHUFFLE(2,2,2,2)));

  return 1;
}

#ifdef DEBUG
static int printf_m128(__m128 v) {
    /* Print the contents of a SSE float4 vector to stdout (debugging) */
  float* p = (float*)(&v);
  printf("%f %f %f %f\n", p[0], p[1], p[2], p[3]);
  return 1;
}
#endif

#endif /* _SSE_TOOLS_H */
