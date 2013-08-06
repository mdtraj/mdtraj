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

#include <stdio.h>
#include <math.h>
#include <pmmintrin.h>
#ifdef HAVE_SSE4
#include <smmintrin.h>
#endif

/****************************************************************************/
/* Utilities                                                                */
/****************************************************************************/

static inline __m128 load_float3(const float* value) {
  // Load (x,y,z) into a SSE register, leaving the last entry
  // set to zero.
  __m128 x = _mm_load_ss(&value[0]);
  __m128 y = _mm_load_ss(&value[1]);
  __m128 z = _mm_load_ss(&value[2]);
  __m128 xy = _mm_movelh_ps(x, y);
  return _mm_shuffle_ps(xy, z, _MM_SHUFFLE(2, 0, 2, 0));
}

static int printf_m128(__m128 v) {
  float* p = (float*)(&v);
  printf("%f %f %f %f\n", p[0], p[1], p[2], p[3]);
  return 1;
}


static int inverse33(const float M[9], __m128 cols[3]) {
  /* Compute the inverse of a 3x3 matrix, storing the columns of the
   * result into three __m128 SSE registers
   */
  double det = M[0] * (M[4] * M[8] - M[5] * M[7])
             + M[3] * (M[7] * M[2] - M[8] * M[1])
             + M[6] * (M[1] * M[5] - M[2] * M[4]);
  __m128 inverse_det = _mm_set1_ps((float) (1.0 / det));
  // cols[0] = _mm_mul_ps(inverse_det, _mm_setr_ps(
  //       M[4]*M[8] - M[7]*M[5], -(M[1]*M[8] - M[2]*M[7]),
  //       M[1]*M[5] - M[2]*M[4], 0.0f));
  // cols[1] = _mm_mul_ps(inverse_det, _mm_setr_ps(
  //       -(M[3]*M[8] - M[5]*M[6]),  M[0]*M[8] - M[2]*M[6],
  //       -(M[0]*M[5] - M[3]*M[2]) , 0.0f));
  // cols[2] = _mm_mul_ps(inverse_det, _mm_setr_ps(
  //       M[3]*M[7] - M[6]*M[4] , -(M[0]*M[7] - M[6]*M[1]),
  //       M[0]*M[4] - M[3]*M[1] , 0.0f));
  //       
  //       
  cols[0] = _mm_mul_ps(inverse_det, _mm_setr_ps(
       M[4]*M[8] - M[7]*M[5], -(M[3]*M[8] - M[5]*M[6]), 
       M[3]*M[7] - M[6]*M[4], 0.0f));
  cols[1] = _mm_mul_ps(inverse_det, _mm_setr_ps(
     -(M[1]*M[8] - M[2]*M[7]), M[0]*M[8] - M[2]*M[6],
     -(M[0]*M[7] - M[6]*M[1]), 0.0f));
   cols[2] = _mm_mul_ps(inverse_det, _mm_setr_ps(
      M[1]*M[5] - M[2]*M[4], -(M[0]*M[5] - M[3]*M[2]),
      M[0]*M[4] - M[3]*M[1] , 0.0f));
 
  return 1;
}

/****************************************************************************/
/* Distances Kernels                                                        */
/****************************************************************************/

int dist(const float* xyz, const int* pairs, float* distance_out,
         float* displacement_out, const int n_frames, const int n_atoms,
         const int n_pairs) {
  /* Compute the distance/displacement  between pairs of atoms in every frame
     of xyz.

     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     pairs : array, shape=(n_pairs, 2)
         The specific pairs of atoms whose distance you want to compute. A 2d
         array of pairs, in C order.
     distance_out : array, shape=(n_frames, n_pairs), optional
         Array where the distances between pairs will be stored, in contiguous
         C order. If NULL is passed in, this return value will not be saved
     displacement_out : array, shaoe=(n_frames, n_pairs, 3), optional
         An optional return value: if you'd also like to save the displacement
         vectors between the pairs, you can pass a pointer here. If
         displacement_out is NULL, then this variable will not be saved back
         to memory.

     All of the arrays are assumed to be contiguous. This code will
     segfault if they're not.
   */

  int i, j;
  int store_displacement = displacement_out == NULL ? 0 : 1;
  int store_distance = distance_out == NULL ? 0 : 1;
  __m128 x1, x2, r12, r12_2, s;

  for (i = 0; i < n_frames; i++) {
    for (j = 0; j < n_pairs; j++) {
      // Load the two vectors whos distance we want to compute
      // x1 = xyz[i, pairs[j,0], 0:3]
      // x2 = xyz[i, pairs[j,1], 0:3]
      x1 = load_float3(xyz + 3*pairs[2*j + 0]);
      x2 = load_float3(xyz + 3*pairs[2*j + 1]);

      // r12 = x2 - x1
      r12 = _mm_sub_ps(x2, x1);
      // r12_2 = r12*r12
      r12_2 = _mm_mul_ps(r12, r12);
      
      if (store_displacement) {
        // store the two lower entries (x,y) in memory
        _mm_storel_pi((__m64*)(displacement_out), r12);
        displacement_out += 2;
        // swap high-low and then store the z entry in the memory
        _mm_store_ss(displacement_out++, _mm_movehl_ps(r12, r12));
      }
      if (store_distance) { 
        // horizontal add the components of d2 with
        // two instructions. note: it's critical
        // here that the last entry of x1 and x2 was 0
        // so that d2.w = 0
        s = _mm_hadd_ps(r12_2, r12_2);
        s = _mm_hadd_ps(s, s);
        // sqrt our final answer
        s = _mm_sqrt_ps(s);

        // s now contains our answer in all four elements, because
        // of the way the hadd works. we only want to store one
        // element.
        _mm_store_ss(distance_out++, s);
      }
    }

    // advance to the next frame
    xyz += n_atoms*3;
  }

  return 1;
}


int dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
             float* distance_out, float* displacement_out,
             const int n_frames, const int n_atoms, const int n_pairs) {
  /* Compute the distance/displacement between pairs of atoms in every frame
     of xyz following the minimum image convention in periodic boundary
     conditions.
    
    The computation follows scheme B.9 in Tukerman, M. "Statistical
    Mechanics: Theory and Molecular Simulation", 2010.

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
     distance_out : array, shape=(n_frames, n_pairs)
         Array where the distances between pairs will be stored, in contiguous
         C order.
     displacement_out : array, shaoe=(n_frames, n_pairs, 3), optional
         An optional return value: if you'd also like to save the displacement
         vectors between the pairs, you can pass a pointer here. If
         displacement_out is NULL, then this variable will not be saved back
         to memory.
         
     All of the arrays are assumed to be contiguous. This code will
     segfault if they're not.
  */
#ifndef HAVE_SSE4
   _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
   int rounding_mode = _MM_GET_ROUNDING_MODE()
#endif

  int i, j;
  int store_displacement = displacement_out == NULL ? 0 : 1;
  int store_distance = distance_out == NULL ? 0 : 1;

  __m128 r1, r2, s1, s2, s12, r12, s, r12_2;
  __m128 hinv[3];
  __m128 h[3];

  for (i = 0; i < n_frames; i++) {
    // Store the columns of the box matrix in three float4s. This format
    // is fast for matrix * vector product. See, for example, this S.O. question:
    // http://stackoverflow.com/questions/14967969/efficient-4x4-matrix-vector-multiplication-with-sse-horizontal-add-and-dot-prod
    h[0] = _mm_setr_ps(box_matrix[0], box_matrix[3], box_matrix[6], 0.0f);
    h[1] = _mm_setr_ps(box_matrix[1], box_matrix[4], box_matrix[7], 0.0f);
    h[2] = _mm_setr_ps(box_matrix[2], box_matrix[5], box_matrix[8], 0.0f);
    // Calculate the inverse of the box matrix, and also store it in the same
    // format.
    inverse33(box_matrix, hinv);
    
    for (j = 0; j < n_pairs; j++) {
      // Load the two vectors whos distance we want to compute
      r1 = load_float3(xyz + 3*pairs[2*j + 0]);
      r2 = load_float3(xyz + 3*pairs[2*j + 1]);
      r12 =  _mm_sub_ps(r2, r1);

      // s12 = INVERSE(H) * r12
      s12 = _mm_add_ps(_mm_add_ps(
         _mm_mul_ps(hinv[0], _mm_shuffle_ps(r12, r12, _MM_SHUFFLE(0,0,0,0))),
         _mm_mul_ps(hinv[1], _mm_shuffle_ps(r12, r12, _MM_SHUFFLE(1,1,1,1)))),
         _mm_mul_ps(hinv[2], _mm_shuffle_ps(r12, r12, _MM_SHUFFLE(2,2,2,2))));

      // s12 = s12 - NEAREST_INTEGER(s12)
#ifdef HAVE_SSE4
      s12 = _mm_sub_ps(s12, _mm_round_ps(s12, _MM_FROUND_TO_NEAREST_INT));
#else
      s12 = _mm_sub_ps(s12, _mm_cvtepi32_ps(_mm_cvtps_epi32(s12)));
#endif

      r12 = _mm_add_ps(_mm_add_ps(
          _mm_mul_ps(h[0], _mm_shuffle_ps(s12, s12, _MM_SHUFFLE(0,0,0,0))),
          _mm_mul_ps(h[1], _mm_shuffle_ps(s12, s12, _MM_SHUFFLE(1,1,1,1)))),
          _mm_mul_ps(h[2], _mm_shuffle_ps(s12, s12, _MM_SHUFFLE(2,2,2,2))));

      if (store_displacement) {
        // store the two lower entries (x,y) in memory
        _mm_storel_pi((__m64*)(displacement_out), r12);
        displacement_out += 2;
        // swap high-low and then store the z entry in the memory
        _mm_store_ss(displacement_out++, _mm_movehl_ps(r12, r12));
      }
      if (store_distance) {
        // out = sqrt(sum(r12**2))
        r12_2 = _mm_mul_ps(r12, r12);
        s = _mm_hadd_ps(r12_2, r12_2);
        s = _mm_hadd_ps(s, s);
        s = _mm_sqrt_ps(s);
        _mm_store_ss(distance_out++, s);
      }
    }
    // advance to the next frame
    xyz += n_atoms*3;
    box_matrix += 9;
  }

#ifndef HAVE_SSE4
   _MM_SET_ROUNDING_MODE(rounding_mode);
#endif
  return 1;
}

/****************************************************************************/
/* Angle Kernels                                                            */
/****************************************************************************/

int angle(const float* xyz, const int* triplets, float* out,
          const int n_frames, const int n_atoms, const int n_pairs) {
  /* Compute the angle between tripples of atoms in every frame
     of xyz.
    
     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     triplets : array, shape=(n_pairs, 3)
         The specific tripple of atoms whose angle you want to compute. The
         angle computed will be centered around the middle element (i.e aABC).
         A 2d array of indices, in C order.
     out : array, shape=(n_frames, n_pairs)
         Array where the angles will be stored, in contiguous C order.

     All of the arrays are assumed to be contiguous. This code will
     segfault if they're not.
     
     Some of the SSE code is adapted from "FastC++: Coding Cpp Efficiently Source Examples"
     Copyright (C) 2011-2013 Matthias Straka, and licensed under the GNU LGPL3.0.
     http://fastcpp.blogspot.com/2011/12/simple-vector3-class-with-sse-support.html
     
     Thanks Matthias!
  */
  
  int i, j;
  __m128 r_m, r_n, r_o, u_prime, u, v_prime, v;
  
  for (i = 0; i < n_frames; i++) {
    for (j = 0; j < n_pairs; j++) {
      r_m = load_float3(xyz + 3*triplets[3*j + 0]);
      r_n = load_float3(xyz + 3*triplets[3*j + 1]);
      r_o = load_float3(xyz + 3*triplets[3*j + 2]);
    
      u_prime = _mm_sub_ps(r_m, r_o);
      v_prime = _mm_sub_ps(r_n, r_o);

      // normalize the vectors u_prime and v_prime
      u = _mm_mul_ps(u_prime, _mm_rsqrt_ps(_mm_dp_ps(u_prime, u_prime, 0x7F)));
      v = _mm_mul_ps(v_prime, _mm_rsqrt_ps(_mm_dp_ps(v_prime, v_prime, 0x7F)));
      
      // compute the arccos of the dot product, and store the result.
      *(out++) = acos(_mm_cvtss_f32(_mm_dp_ps(u_prime, v_prime, 0x71)));
    }
    // advance to the next frame
    xyz += n_atoms*3;
  }
}