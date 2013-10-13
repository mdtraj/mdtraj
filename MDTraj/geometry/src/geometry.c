/* This file is part of MDTraj
 *
 * Copyright 2013 Stanford University
 *
 * MDTraj is free software; you can redistribute it and/or modify
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

#include <stdlib.h>
#include <math.h>
#include <pmmintrin.h>
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#include "ssetools.h"
#include "msvccompat.h"

#ifdef _MSC_VER
 #include "float.h"
 #define isnan(x) (_isnan(x))
#endif


/****************************************************************************/
/* Utilities                                                                */
/****************************************************************************/

static int inverse33(const float M[9], _m128 cols[3]) {
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

INLINE __m128 cross(const __m128 a, const __m128 b) {
   return _mm_sub_ps(
    _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 1, 0, 2))),
    _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1)))
   );
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
#ifndef __SSE4_1__
   _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
   int rounding_mode = _MM_GET_ROUNDING_MODE();
#endif

  int i, j;
  int store_displacement = displacement_out == NULL ? 0 : 1;
  int store_distance = distance_out == NULL ? 0 : 1;

  __m128 r1, r2, s12, r12, s, r12_2;
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
    inverse33(box_matrix, h);
    
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
#ifdef __SSE4_1__
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

#ifndef __SSE4_1__
   _MM_SET_ROUNDING_MODE(rounding_mode);
#endif
  return 1;
}


/****************************************************************************/
/* Angle Kernels                                                            */
/****************************************************************************/

int angle(const float* xyz, const int* triplets, float* out,
          const int n_frames, const int n_atoms, const int n_angles) {
  /* Compute the angle between tripples of atoms in every frame
     of xyz.

     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     triplets : array, shape=(n_angles, 3)
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
    for (j = 0; j < n_angles; j++) {
      r_m = load_float3(xyz + 3*triplets[3*j + 0]);
      r_o = load_float3(xyz + 3*triplets[3*j + 1]);
      r_n = load_float3(xyz + 3*triplets[3*j + 2]);

      u_prime = _mm_sub_ps(r_m, r_o);
      v_prime = _mm_sub_ps(r_n, r_o);

      // normalize the vectors u_prime and v_prime
      u = _mm_mul_ps(u_prime, _mm_rsqrt_ps(_mm_dp_ps(u_prime, u_prime, 0x7F)));
      v = _mm_mul_ps(v_prime, _mm_rsqrt_ps(_mm_dp_ps(v_prime, v_prime, 0x7F)));

      // compute the arccos of the dot product, and store the result.
      *(out++) = acos(_mm_cvtss_f32(_mm_dp_ps(u, v, 0x71)));
    }
    // advance to the next frame
    xyz += n_atoms*3;
  }

  return 1;
}


/****************************************************************************/
/* Dihedral Kernels                                                         */
/****************************************************************************/

int dihedral(const float* xyz, const int* quartets, float* out,
          const int n_frames, const int n_atoms, const int n_quartets) {
  /* Compute the angle between sets of four atoms in every frame
     of xyz.

     Parameters
     ----------
     xyz : array, shape=(n_frames, n_atoms, 3)
         Cartesian coordinates of the atoms in every frame, in contiguous C order.
     quartets : array, shape=(n_quartets, 3)
         The specific quartet of atoms whose angle you want to compute. The
         angle computed will be the torsion around the bound between the
         middle two elements (i.e aABCD). A 2d array of indices, in C order.
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
  __m128 x0, x1, x2, x3, b1, b2, b3, c1, c2, p1, p2;

  for (i = 0; i < n_frames; i++) {
    for (j = 0; j < n_quartets; j++) {
      x0 = load_float3(xyz + 3*quartets[4*j + 0]);
      x1 = load_float3(xyz + 3*quartets[4*j + 1]);
      x2 = load_float3(xyz + 3*quartets[4*j + 2]);
      x3 = load_float3(xyz + 3*quartets[4*j + 3]);

      b1 = _mm_sub_ps(x1, x0);
      b2 = _mm_sub_ps(x2, x1);
      b3 = _mm_sub_ps(x3, x2);

      c1 = cross(b2, b3);
      c2 = cross(b1, b2);

      p1 = _mm_mul_ps(_mm_dp_ps(b1, c1, 0x71), _mm_sqrt_ps(_mm_dp_ps(b2, b2, 0x71)));
      p2 = _mm_dp_ps(c1, c2, 0x71);

      *(out++) = atan2(_mm_cvtss_f32(p1), _mm_cvtss_f32(p2));
    };
    xyz += n_atoms*3;
  }
  return 1;
}


/****************************************************************************/
/* HBond Kernels                                                            */
/****************************************************************************/

static float ks_donor_acceptor(const float* xyz, const float* hcoords,
			       const int* nco_indices, int donor, int acceptor)
{
  /* Conpute the Kabsch-Sander hydrogen bond energy between two residues
     in a single conformation.

     Parameters
     ----------
     xyz : array, shape=(n_atoms, 3)
         All of the atoms in this frame
     nhco0 : array, shape=(4,)
         The indices of the backbone N, H, C, and O atoms in one residue.
     nhco1 : array, shape=(4,)
         The indices of the backbone N, H, C, and O atoms in the other residue.
     donor : int
         Boolean flag. If 0, then nhco0 is the hydrogen bond proton donor (i.e. we
         look at its N and H). If 1, then nhco1 is the hydrogen bond proton donor.

     Returns
     -------
     energy : float
         The KS backbone hydrogen bond energy, in kcal/mol. A number under -0.5
         is considered significant.
  */
  float energy;
  __m128 r_n, r_h, r_c, r_o, r_ho, r_nc, r_hc, r_no, d2_honchcno;
  __m128 coupling;

  // 332 (kcal*A/mol) * 0.42 * 0.2 * (1nm / 10 A)
  coupling = _mm_setr_ps(-2.7888, -2.7888, 2.7888, 2.7888);
  r_n = load_float3(xyz + 3*nco_indices[3*donor]);
  r_h = load_float3(hcoords + 3*donor);
  r_c = load_float3(xyz + 3*nco_indices[3*acceptor + 1]);
  r_o = load_float3(xyz + 3*nco_indices[3*acceptor + 2]);

  //printf("Donor Index %d\n", donor);
  //printf("Acceptor Index %d\n", acceptor);
  /*printf("N index %d\n", 3*nco_indices[3*donor + 0]);
  printf("C index %d\n", 3*nco_indices[3*acceptor + 1]);
  printf("O index %d\n", 3*nco_indices[3*acceptor + 2]);
  printf("\nrN ");
  printf_m128(r_n);
  printf("rH ");
  printf_m128(r_h);
  printf("rC ");
  printf_m128(r_c);
  printf("rO ");
  printf_m128(r_o);*/

  r_ho = _mm_sub_ps(r_h, r_o);
  r_hc = _mm_sub_ps(r_h, r_c);
  r_nc = _mm_sub_ps(r_n, r_c);
  r_no = _mm_sub_ps(r_n, r_o);

  // compute all four dot products (each of the squared distances), and then
  // pack them into a single float4 using three shuffles.
  d2_honchcno = _mm_shuffle_ps(_mm_shuffle_ps(_mm_dp_ps(r_ho, r_ho, 0xF3), _mm_dp_ps(r_nc, r_nc, 0xF3), _MM_SHUFFLE(0,1,0,1)),
                               _mm_shuffle_ps(_mm_dp_ps(r_hc, r_hc, 0xF3), _mm_dp_ps(r_no, r_no, 0xF3), _MM_SHUFFLE(0,1,0,1)),
                               _MM_SHUFFLE(2,0,2,0));

  energy = _mm_cvtss_f32(_mm_dp_ps(coupling, _mm_rsqrt_ps(d2_honchcno), 0xFF));
  //printf("Energy: %f\n\n", energy);
  return (energy < -9.9f ? -9.9f : energy);
}


static int ks_assign_hydrogens(const float* xyz, const int* nco_indices, const int n_residues, float *hcoords)
/* Assign hydrogen atom coordinates
 */
{
  int ri, pc_index, po_index;
  __m128 pc, po, r_co, r_h, r_n, norm_r_co;
  __m128 tenth = _mm_set1_ps(0.1f);

  r_n = load_float3(xyz + 3*nco_indices[0]);
  store_float3(hcoords, r_n);
  hcoords += 3;

  for (ri = 1; ri < n_residues; ri++) {
    pc_index = nco_indices[3*(ri-1) + 1];
    po_index = nco_indices[3*(ri-1) + 2];

    pc = load_float3(xyz + 3*pc_index);
    po = load_float3(xyz + 3*po_index);
    r_co = _mm_sub_ps(pc, po);
    r_n = load_float3(xyz + 3*nco_indices[3*ri + 0]);
    norm_r_co = _mm_mul_ps(r_co, _mm_rsqrt_ps(_mm_dp_ps(r_co, r_co, 0xFF)));
    r_h = _mm_add_ps(r_n, _mm_mul_ps(tenth, norm_r_co));
    store_float3(hcoords, r_h);
    hcoords += 3;
  }

  return 1;
}


static INLINE void store_energies(int* hbonds, float* henergies, int donor,
                             int acceptor, float e) {
  /* Store a computed hbond energy and the appropriate residue indices
     in the output arrays. This function is called twice by kabsch_sander,
     so it seemed appropriate to factor it out.
  */

  float existing_e0 = henergies[2*acceptor + 0];
  float existing_e1 = henergies[2*acceptor + 1];

  if (isnan(existing_e0) || e < existing_e0) {
    //copy over any info in #0 hbond to #1
    hbonds[2*acceptor + 1] = hbonds[acceptor*2 + 0];
    henergies[2*acceptor + 1] = existing_e0;
    hbonds[2*acceptor + 0] = donor;
    henergies[2*acceptor + 0] = e;
    // printf("hbond being stored from donor=%d to acceptor=%d\n", donor, acceptor);
  } else if (isnan(existing_e1) || e < henergies[2*acceptor + 1]) {
    hbonds[2*acceptor + 1] = donor;
    henergies[2*acceptor + 1] = e;
    // printf("hbond being stored from donor=%d to acceptor=%d\n", donor, acceptor);
  }
}

int kabsch_sander(const float* xyz, const int* nco_indices, const int* ca_indices,
                  const int n_frames, const int n_atoms, const int n_residues,
                  int* hbonds, float* henergies) {
  /* Find all of backbone hydrogen bonds between residues in each frame of a
     trajectory.

    Parameters
    ----------
    xyz : array, shape=(n_frames, n_atoms, 3)
        The cartesian coordinates of all of the atoms in each frame.
    nco_indices : array, shape=(n_residues, 3)
        The indices of the backbone N, C, and O atoms for each residue.
    ca_indices : array, shape=(n_residues,)
        The index of the CA atom of each residue. If a residue does not contain
        a CA atom, or you want to skip the residue for another reason, the
	value should be -1

    Returns
    -------
    hbonds : array, shape=(n_frames, n_residues, 2)
        This is a little tricky, so bear with me. This array gives the indices
        of the residues that each backbone hbond *acceptor* is engaged in an hbond
        with. For instance, the equality `bonds[i, j, 0] == k` is interpreted as
        "in frame i, residue j is accepting its first hydrogen bond from residue
        k". `bonds[i, j, 1] == k` means that residue j is accepting its second
        hydrogen bond from residue k. A negative value indicates that no such
        hbond exists.
    henergies : array, shape=(n_frames, n_residues, 2)
        The semantics of this array run parallel to the hbonds array, but
        instead of giving the identity of the interaction partner, it gives
        the energy of the hbond. Only hbonds with energy below -0.5 kcal/mol
        are recorded.
  */

  int i, ri, rj;
  static float HBOND_ENERGY_CUTOFF = -0.5;
  __m128 ri_ca, rj_ca, r12;
  __m128 MINIMAL_CA_DISTANCE2 = _mm_set1_ps(0.81);
  float* hcoords = (float*) malloc(n_residues*3 * sizeof(float));
  if (hcoords == NULL) {
    fprintf(stderr, "Memory Error\n");
    exit(1);
  }

  for (i = 0; i < n_frames; i++) {
    ks_assign_hydrogens(xyz, nco_indices, n_residues, hcoords);

    for (ri = 0; ri < n_residues; ri++) {
      // -1 is used to indicate that this residue lacks a this atom type
      // so just skip it
      if (ca_indices[ri] == -1) continue;
      ri_ca = load_float3(xyz + 3*ca_indices[ri]);

      for (rj = ri + 1; rj < n_residues; rj++) {
        if (ca_indices[rj] == -1) continue;
        rj_ca = load_float3(xyz + 3*ca_indices[rj]);

        // check the ca distance before proceding
        r12 = _mm_sub_ps(ri_ca, rj_ca);
        if(_mm_extract_epi16(CAST__M128I(_mm_cmplt_ps(_mm_dp_ps(r12, r12, 0x7F), MINIMAL_CA_DISTANCE2)), 0)) {
          float e = ks_donor_acceptor(xyz, hcoords, nco_indices, ri, rj);
          if (e < HBOND_ENERGY_CUTOFF)
            // hbond from donor=ri to acceptor=rj
            store_energies(hbonds, henergies, ri, rj, e);

          if (rj != ri + 1) {
	    float e = ks_donor_acceptor(xyz, hcoords, nco_indices, rj, ri);
            if (e < HBOND_ENERGY_CUTOFF)
              // hbond from donor=rj to acceptor=ri
              store_energies(hbonds, henergies, rj, ri, e);
          }
        }
      }
    }
    xyz += n_atoms*3; // advance to the next frame
    hbonds += n_residues*2;
    henergies += n_residues*2;
  }
  free(hcoords);
  return 1;
}
