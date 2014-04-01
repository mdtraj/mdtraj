#ifndef _SSE_TOOLS_H_
#define _SSE_TOOLS_H_

#include "msvccompat.h"
#include <stdio.h>
#include <pmmintrin.h>
#include <smmintrin.h>

/****************************************************************************/
/* Utilities                                                                */
/****************************************************************************/

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

static int printf_m128(__m128 v) {
    /* Print the contents of a SSE float4 vector to stdout (debugging) */
  float* p = (float*)(&v);
  printf("%f %f %f %f\n", p[0], p[1], p[2], p[3]);
  return 1;
}

#endif /* _SSE_TOOLS_H */
