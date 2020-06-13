#ifndef OPENMM_VECTORIZE_SSE_H_
#define OPENMM_VECTORIZE_SSE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

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


static inline __m128 _mm_hsum_ps(__m128 v) {
    v = _mm_hadd_ps(v, v);
    v = _mm_hadd_ps(v, v);
    return v;
}


static inline __m128 _mm_round_ps2(const __m128 a) {
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

inline __m128 _mm_floor_ps2(const __m128& x) {
    /* http://dss.stephanierct.com/DevBlog/?p=8 */
    __m128i v0 = _mm_setzero_si128();
    __m128i v1 = _mm_cmpeq_epi32(v0,v0);
    __m128i ji = _mm_srli_epi32( v1, 25);
    __m128 j = _mm_castsi128_ps(_mm_slli_epi32( ji, 23)); //create vector 1.0f
    __m128i i = _mm_cvttps_epi32(x);
    __m128 fi = _mm_cvtepi32_ps(i);
    __m128 igx = _mm_cmpgt_ps(fi, x);
    j = _mm_and_ps(igx, j);
    return _mm_sub_ps(fi, j);
}

static inline __m128 _mm_dp_ps2( __m128 a, __m128 b, const int mask ) {
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

// This file defines classes and functions to simplify vectorizing code with SSE.

class ivec4;

/**
 * A four element vector of floats.
 */
class fvec4 {
public:
    __m128 val;

    fvec4() {}
    fvec4(float v) : val(_mm_set1_ps(v)) {}
    fvec4(float v1, float v2, float v3, float v4) : val(_mm_set_ps(v4, v3, v2, v1)) {}
    fvec4(__m128 v) : val(v) {}

    /*  Note on the initialization of fvec4 with a pointer
        --------------------------------------------------------
        The constructor below loads a four-element array. You need to be careful about
        checking that all four elements lead to valid memory, especially when working
        with 3-vectors (xyz). If the fourth element is a NaN or unallocated, you'll get
        incorrect results or crashes.
    */
    fvec4(const float* v) : val(_mm_loadu_ps(v)) {}

    operator __m128() const {
        return val;
    }
    float operator[](int i) const {
        float result[4];
        store(result);
        return result[i];
    }
    void store(float* v) const {
        _mm_storeu_ps(v, val);
    }
    fvec4 operator+(const fvec4& other) const {
        return _mm_add_ps(val, other);
    }
    fvec4 operator-(const fvec4& other) const {
        return _mm_sub_ps(val, other);
    }
    fvec4 operator*(const fvec4& other) const {
        return _mm_mul_ps(val, other);
    }
    fvec4 operator/(const fvec4& other) const {
        return _mm_div_ps(val, other);
    }
    void operator+=(const fvec4& other) {
        val = _mm_add_ps(val, other);
    }
    void operator-=(const fvec4& other) {
        val = _mm_sub_ps(val, other);
    }
    void operator*=(const fvec4& other) {
        val = _mm_mul_ps(val, other);
    }
    void operator/=(const fvec4& other) {
        val = _mm_div_ps(val, other);
    }
    fvec4 operator-() const {
        return _mm_sub_ps(_mm_set1_ps(0.0f), val);
    }
    fvec4 operator&(const fvec4& other) const {
        return _mm_and_ps(val, other);
    }
    fvec4 operator|(const fvec4& other) const {
        return _mm_or_ps(val, other);
    }
    fvec4 operator==(const fvec4& other) const {
        return _mm_cmpeq_ps(val, other);
    }
    fvec4 operator!=(const fvec4& other) const {
        return _mm_cmpneq_ps(val, other);
    }
    fvec4 operator>(const fvec4& other) const {
        return _mm_cmpgt_ps(val, other);
    }
    fvec4 operator<(const fvec4& other) const {
        return _mm_cmplt_ps(val, other);
    }
    fvec4 operator>=(const fvec4& other) const {
        return _mm_cmpge_ps(val, other);
    }
    fvec4 operator<=(const fvec4& other) const {
        return _mm_cmple_ps(val, other);
    }
    operator ivec4() const;
};

/**
 * A four element vector of ints.
 */
class ivec4 {
public:
    __m128i val;

    ivec4() {}
    ivec4(int v) : val(_mm_set1_epi32(v)) {}
    ivec4(int v1, int v2, int v3, int v4) : val(_mm_set_epi32(v4, v3, v2, v1)) {}
    ivec4(__m128i v) : val(v) {}
    ivec4(const int* v) : val(_mm_loadu_si128((const __m128i*) v)) {}
    operator __m128i() const {
        return val;
    }
    int operator[](int i) const {
        int result[4];
        store(result);
        return result[i];
    }
    void store(int* v) const {
        _mm_storeu_si128((__m128i*) v, val);
    }
    ivec4 operator+(const ivec4& other) const {
        return _mm_add_epi32(val, other);
    }
    ivec4 operator-(const ivec4& other) const {
        return _mm_sub_epi32(val, other);
    }
    /*
    ivec4 operator*(const ivec4& other) const {
        return _mm_mullo_epi32(val, other);
    }
    */
    void operator+=(const ivec4& other) {
        val = _mm_add_epi32(val, other);
    }
    void operator-=(const ivec4& other) {
        val = _mm_sub_epi32(val, other);
    }
    /*
    void operator*=(const ivec4& other) {
        val = _mm_mullo_epi32(val, other);
    }
    */
    ivec4 operator-() const {
        return _mm_sub_epi32(_mm_set1_epi32(0), val);
    }
    ivec4 operator&(const ivec4& other) const {
        return _mm_and_si128(val, other);
    }
    ivec4 operator|(const ivec4& other) const {
        return _mm_or_si128(val, other);
    }
    ivec4 operator==(const ivec4& other) const {
        return _mm_cmpeq_epi32(val, other);
    }
    ivec4 operator!=(const ivec4& other) const {
        return _mm_xor_si128(*this==other, _mm_set1_epi32(0xFFFFFFFF));
    }
    ivec4 operator>(const ivec4& other) const {
        return _mm_cmpgt_epi32(val, other);
    }
    ivec4 operator<(const ivec4& other) const {
        return _mm_cmplt_epi32(val, other);
    }
    ivec4 operator>=(const ivec4& other) const {
        return _mm_xor_si128(_mm_cmplt_epi32(val, other), _mm_set1_epi32(0xFFFFFFFF));
    }
    ivec4 operator<=(const ivec4& other) const {
        return _mm_xor_si128(_mm_cmpgt_epi32(val, other), _mm_set1_epi32(0xFFFFFFFF));
    }
    operator fvec4() const;
};

// Conversion operators.

inline fvec4::operator ivec4() const {
    return _mm_cvttps_epi32(val);
}

inline ivec4::operator fvec4() const {
    return _mm_cvtepi32_ps(val);
}

// Functions that operate on fvec4s.

static inline fvec4 round(const fvec4& v) {
    return fvec4(_mm_round_ps2(v.val));
}

static inline fvec4 floor(const fvec4& v) {
    return fvec4(_mm_floor_ps2(v.val));
}

static inline fvec4 min(const fvec4& v1, const fvec4& v2) {
    return fvec4(_mm_min_ps(v1.val, v2.val));
}

static inline fvec4 max(const fvec4& v1, const fvec4& v2) {
    return fvec4(_mm_max_ps(v1.val, v2.val));
}

static inline fvec4 abs(const fvec4& v) {
    static const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
    return fvec4(_mm_and_ps(v.val, mask));
}

static inline fvec4 sqrt(const fvec4& v) {
    return fvec4(_mm_sqrt_ps(v.val));
}

static inline fvec4 rsqrt(const fvec4& v) {
    // Initial estimate of rsqrt().

    fvec4 y(_mm_rsqrt_ps(v.val));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = v*0.5f;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}

static inline float dot3(const fvec4& v1, const fvec4& v2) {
    return _mm_cvtss_f32(_mm_dp_ps2(v1, v2, 0x71));
}

static inline float dot4(const fvec4& v1, const fvec4& v2) {
    return _mm_cvtss_f32(_mm_dp_ps2(v1, v2, 0xF1));
}

static inline fvec4 cross(const fvec4& v1, const fvec4& v2) {
    fvec4 temp = fvec4(_mm_mul_ps(v1, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 0, 2, 1)))) -
                 fvec4(_mm_mul_ps(v2, _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1))));
    return _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(3, 0, 2, 1));
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4) {
    _MM_TRANSPOSE4_PS(v1, v2, v3, v4);
}

// Mathematical operators involving a scalar and a vector.

static inline fvec4 operator+(float v1, const fvec4& v2) {
    return fvec4(v1)+v2;
}

static inline fvec4 operator-(float v1, const fvec4& v2) {
    return fvec4(v1)-v2;
}

static inline fvec4 operator*(float v1, const fvec4& v2) {
    return fvec4(v1)*v2;
}

static inline fvec4 operator/(float v1, const fvec4& v2) {
    return fvec4(v1)/v2;
}

static inline fvec4 load3(const float* v) {
    /* Load (x,y,z) into a SSE register, leaving the last entry */
    /* set to zero. */
  __m128 x = _mm_load_ss(&v[0]);
  __m128 y = _mm_load_ss(&v[1]);
  __m128 z = _mm_load_ss(&v[2]);
  __m128 xy = _mm_movelh_ps(x, y);
  return fvec4(_mm_shuffle_ps(xy, z, _MM_SHUFFLE(2, 0, 2, 0)));
}


static inline int store3(const fvec4& v, float* loc) {
    /* Store the low three floats in an SSE register into */
    /* memory, at location loc[0], loc[1], loc[2]. The high */
    /* float is not touched. */
  __m128 val = v.val;
  _mm_store_ss(loc, val);
  _mm_store_ss(loc+1, _mm_shuffle_ps(val, val, _MM_SHUFFLE(1,1,1,1)));
  _mm_store_ss(loc+2, _mm_shuffle_ps(val, val, _MM_SHUFFLE(2,2,2,2)));

  return 1;
}


#endif /*OPENMM_VECTORIZE_SSE_H_*/

