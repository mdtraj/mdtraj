/* This file is part of MDTraj.
 *
 * Copyright 2013 Stanford University
 *
 * MSMBuilder is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <pmmintrin.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <stdio.h>


static inline __m128 load_float3(float* value) {
  // Load (x,y,z) into a SSE register, leaving the last entry
  // set to zero.
  __m128 x = _mm_load_ss(&value[0]);
  __m128 y = _mm_load_ss(&value[1]);
  __m128 z = _mm_load_ss(&value[2]);
  __m128 xy = _mm_movelh_ps(x, y);
  return _mm_shuffle_ps(xy, z, _MM_SHUFFLE(2, 0, 2, 0));
}

int dist(float* xyz, int* pairs, float* out,
         int n_frames, int n_atoms, int n_pairs) {
  /* Compute the distance between pairs of atoms in every frame
     of xyz.

     Uses SSE3 horizontal add for the reduction (compile with -msse3)

     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     pairs : array, shape=(n_pairs, 2)
         The specific pairs of atoms whose distance you want to compute. A 2d
         array of pairs, in C order.
     out : array, shape=(n_frames, n_pairs)
         Array where the output will be stored, in contiguous C order.

     All of the arrays are assumed to be contiguous. This code will
     segfault if they're not.
   */

  int i, j;
  __m128 x1, x2, r, r2, s;

  for (i = 0; i < n_frames; i++) {
    for (j = 0; j < n_pairs; j++) {
      // Load the two vectors whos distance we want to compute
      // x1 = xyz[i, pairs[j,0], 0:3]
      // x2 = xyz[i, pairs[j,1], 0:3]
      x1 = load_float3(xyz + 3*pairs[2*j + 0]);
      x2 = load_float3(xyz + 3*pairs[2*j + 1]);

      // r = x2 - x1
      r = _mm_sub_ps(x2, x1);
      // r2 = r*r
      r2 = _mm_mul_ps(r, r);

      // horizontal add the components of d2 with
      // two instructions. note: it's critical
      // here that the last entry of x1 and x2 was 0
      // so that d2.w = 0
      s = _mm_hadd_ps(r2, r2);
      s = _mm_hadd_ps(s, s);
      // sqrt our final answer
      s = _mm_sqrt_ps(s);

      // s now contains our answer in all four elements, because
      // of the way the hadd works. we only want to store one
      // element.
      _mm_store_ss(out++, s);
    }

    // advance to the next frame
    xyz += n_atoms*3;
  }

  return 1;
}


static int printf_m128(__m128 v) {
  float* p = (float*)(&v);
  printf("%f %f %f %f\n", p[0], p[1], p[2], p[3]);
  return 1;
}

static int inverse33(float M[9], __m128 cols[3]) {
  double det = M[0] * (M[4] * M[8] - M[5] * M[7])
             + M[3] * (M[7] * M[2] - M[8] * M[1])
             + M[6] * (M[1] * M[5] - M[2] * M[4]);
  __m128 inverse_det = _mm_set1_ps((float) (1.0 / det));
  cols[0] = _mm_mul_ps(inverse_det, _mm_setr_ps(
        M[4]*M[8] - M[7]*M[5], -(M[1]*M[8] - M[2]*M[7]),
        M[1]*M[5] - M[2]*M[4], 0.0f));
  cols[1] = _mm_mul_ps(inverse_det, _mm_setr_ps(
        -(M[3]*M[8] - M[5]*M[6]),  M[0]*M[8] - M[2]*M[6],
        -(M[0]*M[5] - M[3]*M[2]) , 0.0f));
  cols[2] = _mm_mul_ps(inverse_det, _mm_setr_ps(
        M[3]*M[7] - M[6]*M[4] , -(M[0]*M[7] - M[6]*M[1]),
        M[0]*M[4] - M[3]*M[1] , 0.0f));
  return 1;
}


int dist_mic(float* xyz, int* pairs, float* box_matrix, float* out,
             int n_frames, int n_atoms, int n_pairs) {
  /* Compute the distance between pairs of atoms in every frame
     of xyz following the minimum image convention in periodic boundary
     conditions.

     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     pairs : array, shape=(n_pairs, 2)
         The specific pairs of atoms whose distance you want to compute. A 2d
         array of pairs, in C order.
     box_matrix : array, shape=(3,3)
          The box matrix for a single frame. All of the frames are assumed to
          use this box vector.
     out : array, shape=(n_frames, n_pairs)
         Array where the output will be stored, in contiguous C order.

     All of the arrays are assumed to be contiguous. This code will
     segfault if they're not.
  */

  // Invert the box matrix, and store each row of the result in a float4
  // with zero padding on the last element
  __m128 hinv[3];
  __m128 h[3];
  h[0] = _mm_setr_ps(box_matrix[0], box_matrix[3], box_matrix[6], 0.0f);
  h[1] = _mm_setr_ps(box_matrix[1], box_matrix[4], box_matrix[7], 0.0f);
  h[2] = _mm_setr_ps(box_matrix[2], box_matrix[5], box_matrix[8], 0.0f);
  inverse33(box_matrix, hinv);

  int i, j;
  __m128 r1, r2, s1, s2, s12, r12, s, r12_2;

  for (i = 0; i < n_frames; i++) {
    for (j = 0; j < n_pairs; j++) {
      // Load the two vectors whos distance we want to compute
      r1 = load_float3(xyz + 3*pairs[2*j + 0]);
      r2 = load_float3(xyz + 3*pairs[2*j + 1]);

      // s1 = INVERSE(H) * r1
      s1 = _mm_add_ps(_mm_add_ps(
         _mm_mul_ps(hinv[0], _mm_shuffle_ps(r1, r1, _MM_SHUFFLE(0,0,0,0))),
         _mm_mul_ps(hinv[1], _mm_shuffle_ps(r1, r1, _MM_SHUFFLE(1,1,1,1)))),
         _mm_mul_ps(hinv[2], _mm_shuffle_ps(r1, r1, _MM_SHUFFLE(2,2,2,2))));
      // s2 = INVERSE(H) * r2
      s2 = _mm_add_ps(_mm_add_ps(
         _mm_mul_ps(hinv[0], _mm_shuffle_ps(r2, r2, _MM_SHUFFLE(0,0,0,0))),
         _mm_mul_ps(hinv[1], _mm_shuffle_ps(r2, r2, _MM_SHUFFLE(1,1,1,1)))),
         _mm_mul_ps(hinv[2], _mm_shuffle_ps(r2, r2, _MM_SHUFFLE(2,2,2,2))));

      // s12 = s12 - NEAREST_INTEGER(s12)
      s12 = _mm_sub_ps(s2, s1);
      s12 = _mm_sub_ps(s12, _mm_round_ps(s12, _MM_FROUND_TO_NEAREST_INT));

      r12 = _mm_add_ps(_mm_add_ps(
          _mm_mul_ps(h[0], _mm_shuffle_ps(s12, s12, _MM_SHUFFLE(0,0,0,0))),
          _mm_mul_ps(h[1], _mm_shuffle_ps(s12, s12, _MM_SHUFFLE(1,1,1,1)))),
          _mm_mul_ps(h[2], _mm_shuffle_ps(s12, s12, _MM_SHUFFLE(2,2,2,2))));

      // out = sqrt(sum(r12**2))
      r12_2 = _mm_mul_ps(r12, r12);
      s = _mm_hadd_ps(r12_2, r12_2);
      s = _mm_hadd_ps(s, s);
      s = _mm_sqrt_ps(s);
      _mm_store_ss(out++, s);
    }
    // advance to the next frame
    xyz += n_atoms*3;
  }


  return 1;
}
