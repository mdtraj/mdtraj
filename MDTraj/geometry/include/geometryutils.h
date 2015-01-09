#ifndef _GEOMETRY_UTILS_H_
#define _GEOMETRY_UTILS_H_

#include "msvccompat.h"
#include "ssetools.h"


/**
 * Compute the cross-product of two 3D vectors stored the first three float
 * elements of `a` and `b`.
 */
static __m128 cross(const __m128 a, const __m128 b)
{
   return _mm_sub_ps(
    _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 1, 0, 2))),
    _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1)))
   );
}


/**
 * Compute the inverse of a 3x3 matrix, storing the columns of the
 * result into three __m128 SSE registers
 */

static int inverse33(const float M[9], __m128* cols_0, __m128* cols_1, __m128* cols_2)
{
  double det = M[0] * (M[4] * M[8] - M[5] * M[7])
             + M[3] * (M[7] * M[2] - M[8] * M[1])
             + M[6] * (M[1] * M[5] - M[2] * M[4]);
  __m128 inverse_det = _mm_set1_ps((float) (1.0 / det));
  *cols_0 = _mm_mul_ps(inverse_det, _mm_setr_ps(
       M[4]*M[8] - M[7]*M[5], -(M[3]*M[8] - M[5]*M[6]),
       M[3]*M[7] - M[6]*M[4], 0.0f));
  *cols_1 = _mm_mul_ps(inverse_det, _mm_setr_ps(
     -(M[1]*M[8] - M[2]*M[7]), M[0]*M[8] - M[2]*M[6],
     -(M[0]*M[7] - M[6]*M[1]), 0.0f));
  *cols_2 = _mm_mul_ps(inverse_det, _mm_setr_ps(
      M[1]*M[5] - M[2]*M[4], -(M[0]*M[5] - M[3]*M[2]),
      M[0]*M[4] - M[3]*M[1] , 0.0f));

  return 1;
}


/**
 * Load a 3x3 matrix in C-order from main memory and store it as well as its
 * matrix inverse in column-major order into SSE __m128 registers.
 *
 * Parameters
 * ----------
 * M : array (input)
 *     3x3 matrix of floats in "c" (row-major) order.
 * h : pointer to array of three __m128s (out)
 *    Each column of the matrix M will be stored in one of the elements of (*h)
 * hinv : pointer to array of three __m128s (out)
 *    Each column of the matrix inverse of M will be stored in one of the
 *    elements of (*hinv)
 *
 * This format is fast for matrix * vector product. See, for example,
 * this S.O. question: http://stackoverflow.com/questions/14967969/efficient-4x4-matrix-vector-multiplication-with-sse-horizontal-add-and-dot-prod
 */
static INLINE void loadBoxMatrix(const float M[9], __m128 (*h)[3], __m128 (*hinv)[3])
{
    (*h)[0] = _mm_setr_ps(M[0], M[3], M[6], 0.0f);
    (*h)[1] = _mm_setr_ps(M[1], M[4], M[7], 0.0f);
    (*h)[2] = _mm_setr_ps(M[2], M[5], M[8], 0.0f);
    /* Calculate the inverse of the box matrix, and also store it in the same */
    /* format. */
    inverse33(M, (*hinv), (*hinv)+1, (*hinv)+2);
}


/**
 * Transform a displacement vector in R^3 in a euclidean space with
 * periodic boundary conditions according to the minimum image convention
 *
 * The math is pretty simple. It can be found in Appendix B, Eq B.9 of [1]
 *
 *     s_ij = H^{-1} * r_ij
 *     s_ij = s_ij - NEAREST_INTEGER(s_ij)
 *     r_ij = H * s_ij
 *
 * where h and h^{-1} are the 3x3 matrix of box vectors and the inverse
 * of that matrix respectively.
 * -------------------------------------------------------------------
 * For this function, r is supplied in the first three elements of float4
 * vector (__m128) `r`. `h` and `hinv` are pointers to arrays of three
 * float4s which specify the matricies H and H^{-1} in column-major order.
 * That is, the first float4 (*h)[0] gives [H[0,0], H[1,0], H[2,0]. 0],
 * (*h)[2] is [H[0,1], H[1,1], H[2,1], 0], etc.
 *
 * .. [1] Tuckerman, Mark. Statistical Mechanics and Molecular Simulations.
 *    Oxford University Press, 2008.
 */
static INLINE __m128 minimum_image(__m128 r, const __m128 (*h)[3], const __m128 (*hinv)[3])
{
    __m128 s;
    s = _mm_add_ps(_mm_add_ps(
       _mm_mul_ps((*hinv)[0], _mm_shuffle_ps(r, r, _MM_SHUFFLE(0,0,0,0))),
       _mm_mul_ps((*hinv)[1], _mm_shuffle_ps(r, r, _MM_SHUFFLE(1,1,1,1)))),
       _mm_mul_ps((*hinv)[2], _mm_shuffle_ps(r, r, _MM_SHUFFLE(2,2,2,2))));

    /* s = s - NEAREST_INTEGER(s) */
    /* s = _mm_sub_ps(s, _mm_round_ps(s, _MM_FROUND_TO_NEAREST_INT)); */
    s = _mm_sub_ps(s, _mm_round_ps2(s));

    r = _mm_add_ps(_mm_add_ps(
        _mm_mul_ps((*h)[0], _mm_shuffle_ps(s, s, _MM_SHUFFLE(0,0,0,0))),
        _mm_mul_ps((*h)[1], _mm_shuffle_ps(s, s, _MM_SHUFFLE(1,1,1,1)))),
        _mm_mul_ps((*h)[2], _mm_shuffle_ps(s, s, _MM_SHUFFLE(2,2,2,2))));
    return r;
}


#endif
