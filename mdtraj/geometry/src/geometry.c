/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2012-2015 Stanford University and the Authors               */
/*                                                                       */
/* Authors: Robert McGibbon                                              */
/* Contributors: Jason Swails                                            */
/*                                                                       */
/* MDTraj is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as        */
/* published by the Free Software Foundation, either version 2.1         */
/* of the License, or (at your option) any later version.                */
/*                                                                       */
/* This library is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU Lesser General Public License for more details.                   */
/*                                                                       */
/* You should have received a copy of the GNU Lesser General Public      */
/* License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.*/
/*=======================================================================*/


#include <stdlib.h>
#include <math.h>
#include "geometry.h"
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define CLIP(X, X_min, X_max) (MIN(MAX(X, X_min), X_max))

#include <pmmintrin.h>
#include "ssetools.h"
#include "msvccompat.h"
#include "geometryutils.h"
#ifdef _MSC_VER
 #include "float.h"
 #define isnan(x) (_isnan(x))
#endif


/****************************************************************************/
/* Distance, Angle and Dihedral kernels                                     */
/****************************************************************************/

/**
 * This is kindof hacky / gross, but I think it's the best way to avoid having
 * lots of copy-paste code. For each of the distance/angle/dihedral kernels, we
 * want to compile two version: one which uses PBCs and the other which does
 * not. Most of the code between these versions is shared, so I've written
 * the parts of the two functions which are different using #ifdefs. So we
 * just include these files _twice_ here, toggling the variable that controls
 * the ifdef.
 *
 * Note that these kernel files are really not capable of being compiled
 * independently -- they're not header files at all -- and they're really just
 * meant to be #included here.
 **/
#undef COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#include "distancekernels.h"
#include "anglekernels.h"
#include "dihedralkernels.h"

#define COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#include "distancekernels.h"
#include "anglekernels.h"
#include "dihedralkernels.h"

/**
 * Distance kernel for general triclinic cells. It is the same as the standard
 * MIC distance code, except that it has to check the distances in all 27
 * surrounding unit cells to pick out the closest one. This is because the MIC
 * code uses scaled, or fractional, coordinates and finds the minimum image that
 * way, which eliminates the anisotropy of the unit cell (i.e., for
 * non-orthorhombic cells, some corners of the unit cell are farther away than
 * others). The only distances that the algorithm in dist_mic are *guaranteed*
 * to get correct are the ones that fall within the largest sphere that can be
 * entirely contained between the planes of the adjacent unit cells (equal to
 * the largest cutoff distance permissible in MD simulations).
 *
 * TODO: Add SSE-vectorization to the nearest neighbor search here
 */
int dist_mic_triclinic(const float* xyz, const int* pairs, const float* box_matrix,
                       float* distance_out, float* displacement_out, int n_frames,
                       int n_atoms, int n_pairs) {
    float *displacements; // We need to get these regardless
    float *distances;     // We need to get these regardless
    int f, i, j, k, n, dist_start, disp_start, dist_idx, disp_idx;
    float min_dist, dist_test;

    float orig_disp[3], min_disp[3], disp_test[3];
    float v1[3], v2[3], bv1[3], bv2[3], bv3[3];

    if (displacement_out == NULL)
        displacements = (float*) malloc(n_pairs*n_frames*3 * sizeof(float));
    else
        displacements = displacement_out;

    if (distance_out == NULL)
        distances = (float*) malloc(n_pairs*n_frames * sizeof(float));
    else
        distances = distance_out;

    dist_mic(xyz, pairs, box_matrix, distances, displacements, n_frames, n_atoms, n_pairs);

    /* Now we have to search the surrounding unit cells  */
    for (f = 0; f < n_frames; f++) {
        dist_start = n_pairs*f;
        disp_start = dist_start*3;

        /* Store the original box vectors, which are columns in the row-major matrix layout */
        bv1[0] = box_matrix[9*f  ]; bv2[0] = box_matrix[9*f+1]; bv3[0] = box_matrix[9*f+2];
        bv1[1] = box_matrix[9*f+3]; bv2[1] = box_matrix[9*f+4]; bv3[1] = box_matrix[9*f+5];
        bv1[2] = box_matrix[9*f+6]; bv2[2] = box_matrix[9*f+7]; bv3[2] = box_matrix[9*f+8];

        for (n = 0; n < n_pairs; n++) {
            dist_idx = dist_start + n;
            disp_idx = disp_start + n*3;
            min_dist = distances[dist_idx]*distances[dist_idx];
            min_disp[0] = displacements[disp_idx  ];
            min_disp[1] = displacements[disp_idx+1];
            min_disp[2] = displacements[disp_idx+2];
            orig_disp[0] = displacements[disp_idx  ];
            orig_disp[1] = displacements[disp_idx+1];
            orig_disp[2] = displacements[disp_idx+2];
            for (i = -1; i < 2; i++) {
                v1[0] = bv1[0]*i;
                v1[1] = bv1[1]*i;
                v1[2] = bv1[2]*i;
                for (j = -1; j < 2; j++) {
                    v2[0] = bv2[0]*j + v1[0];
                    v2[1] = bv2[1]*j + v1[1];
                    v2[2] = bv2[2]*j + v1[2];
                    for (k = -1; k < 2; k++) {
                        disp_test[0] = orig_disp[0] + v2[0] + bv3[0]*k;
                        disp_test[1] = orig_disp[1] + v2[1] + bv3[1]*k;
                        disp_test[2] = orig_disp[2] + v2[2] + bv3[2]*k;
                        dist_test = disp_test[0]*disp_test[0] + disp_test[1]*disp_test[1] + disp_test[2]*disp_test[2];
                        if (dist_test < min_dist) {
                            min_dist = dist_test;
                            min_disp[0] = disp_test[0];
                            min_disp[1] = disp_test[1];
                            min_disp[2] = disp_test[2];
                        }
                    }
                }
            }
            if (distance_out != NULL)
                distances[dist_idx] = sqrtf(min_dist);
            displacements[disp_idx  ] = min_disp[0];
            displacements[disp_idx+1] = min_disp[1];
            displacements[disp_idx+2] = min_disp[2];
        }
    }

    /* If we had to allocate either displacements or distances, deallocate them now */
    if (displacement_out == NULL) free(displacements);
    if (distance_out == NULL) free(distances);

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
  __m128 coupling, recip_sqrt, one;
  one = _mm_set1_ps(1.0);

  /* 332 (kcal*A/mol) * 0.42 * 0.2 * (1nm / 10 A) */
  coupling = _mm_setr_ps(-2.7888, -2.7888, 2.7888, 2.7888);
  r_n = load_float3(xyz + 3*nco_indices[3*donor]);
  r_h = load_float3(hcoords + 3*donor);
  r_c = load_float3(xyz + 3*nco_indices[3*acceptor + 1]);
  r_o = load_float3(xyz + 3*nco_indices[3*acceptor + 2]);

  /*
  printf("Donor Index %d\n", donor);
  printf("Acceptor Index %d\n", acceptor);
  printf("N index %d\n", 3*nco_indices[3*donor + 0]);
  printf("C index %d\n", 3*nco_indices[3*acceptor + 1]);
  printf("O index %d\n", 3*nco_indices[3*acceptor + 2]);
  printf("\nrN ");
  printf_m128(r_n);
  printf("rH ");
  printf_m128(r_h);
  printf("rC ");
  printf_m128(r_c);
  printf("rO ");
  printf_m128(r_o);
  */

  r_ho = _mm_sub_ps(r_h, r_o);
  r_hc = _mm_sub_ps(r_h, r_c);
  r_nc = _mm_sub_ps(r_n, r_c);
  r_no = _mm_sub_ps(r_n, r_o);

  /* compute all four dot products (each of the squared distances), and then */
  /* pack them into a single float4 using three shuffles. */
  d2_honchcno = _mm_shuffle_ps(_mm_shuffle_ps(_mm_dp_ps2(r_ho, r_ho, 0xF3), _mm_dp_ps2(r_nc, r_nc, 0xF3), _MM_SHUFFLE(0,1,0,1)),
                               _mm_shuffle_ps(_mm_dp_ps2(r_hc, r_hc, 0xF3), _mm_dp_ps2(r_no, r_no, 0xF3), _MM_SHUFFLE(0,1,0,1)),
                               _MM_SHUFFLE(2,0,2,0));

  /* rsqrt_ps is really not that accurate... */
  recip_sqrt = _mm_div_ps(one, _mm_sqrt_ps(d2_honchcno));
  energy = _mm_cvtss_f32(_mm_dp_ps2(coupling, recip_sqrt, 0xFF));
  // energy = _mm_cvtss_f32(_mm_dp_ps(coupling, _mm_rsqrt_ps(d2_honchcno), 0xFF));
  return (energy < -9.9f ? -9.9f : energy);
}


static int ks_assign_hydrogens(const float* xyz, const int* nco_indices, const int n_residues, float *hcoords, int* skip)
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
      if (!skip[ri]) {
          pc_index = nco_indices[3*(ri-1) + 1];
          po_index = nco_indices[3*(ri-1) + 2];

          pc = load_float3(xyz + 3*pc_index);
          po = load_float3(xyz + 3*po_index);
          r_co = _mm_sub_ps(pc, po);
          r_n = load_float3(xyz + 3*nco_indices[3*ri + 0]);
          norm_r_co = _mm_mul_ps(r_co, _mm_rsqrt_ps(_mm_dp_ps2(r_co, r_co, 0xFF)));
          r_h = _mm_add_ps(r_n, _mm_mul_ps(tenth, norm_r_co));
          store_float3(hcoords, r_h);
      }
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
  // if (donor == 1 && acceptor == 48)
  // printf("storing donor=1, acceptor=48: energy=%f\n", e);

  float existing_e0 = henergies[2*donor + 0];
  float existing_e1 = henergies[2*donor + 1];

  if (isnan(existing_e0) || e < existing_e0) {
    /* copy over any info in #0 hbond to #1 */
    hbonds[2*donor + 1] = hbonds[donor*2 + 0];
    henergies[2*donor + 1] = existing_e0;
    hbonds[2*donor + 0] = acceptor;
    henergies[2*donor + 0] = e;
    /* printf("hbond being stored from donor=%d to acceptor=%d\n", donor, acceptor); */
  } else if (isnan(existing_e1) || e < henergies[2*donor + 1]) {
    hbonds[2*donor + 1] = acceptor;
    henergies[2*donor + 1] = e;
    /* printf("hbond being stored from donor=%d to acceptor=%d\n", donor, acceptor); */
  }
}

int kabsch_sander(const float* xyz, const int* nco_indices, const int* ca_indices,
                  const int* is_proline, const int n_frames, const int n_atoms,
                  const int n_residues, int* hbonds, float* henergies) {
  /* Find all of backbone hydrogen bonds between residues in each frame of a
     trajectory.

    Parameters
    ----------
    xyz : array, shape=(n_frames, n_atoms, 3)
        The cartesian coordinates of all of the atoms in each frame.
    nco_indices : array, shape=(n_residues, 3)
        The indices of the backbone N, C, and O atoms for each residue.
    ca_indices : array, shape=(n_residues,)
        The index of the CA atom of each residue.
    is_proline : array, shape=(n_residue,)
        If a particular residue does not contain a CA atom, or you want to skip
        the residue for another reason, this value should evaluate to True.

    Returns
    -------
    hbonds : array, shape=(n_frames, n_residues, 2)
        This is a little tricky, so bear with me. This array gives the indices
        of the residues that each backbone hbond *donor* is engaged in an hbond
        with. For instance, the equality `bonds[i, j, 0] == k` is interpreted as
        "in frame i, residue j is donating its first hydrogen bond from residue
        k". `bonds[i, j, 1] == k` means that residue j is donating its second
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
  int* skip = (int*) calloc(n_residues, sizeof(int));
  if (hcoords == NULL || skip == NULL) {
    fprintf(stderr, "Memory Error\n");
    exit(1);
  }
  for (i = 0; i < n_residues; i++)
      if ((nco_indices[i*3] == -1) || (nco_indices[i*3+1] == -1) ||
          (nco_indices[i*3+2] == -1) || ca_indices[i] == -1)
          skip[i] = 1;

  for (i = 0; i < n_frames; i++) {
    ks_assign_hydrogens(xyz, nco_indices, n_residues, hcoords, skip);

    for (ri = 0; ri < n_residues; ri++) {
      if (skip[ri]) continue;
      ri_ca = load_float3(xyz + 3*ca_indices[ri]);

      for (rj = ri + 1; rj < n_residues; rj++) {
        if (skip[rj]) continue;
        rj_ca = load_float3(xyz + 3*ca_indices[rj]);

        /* check the ca distance before proceding */
        r12 = _mm_sub_ps(ri_ca, rj_ca);

        if(_mm_extract_epi16(CAST__M128I(_mm_cmplt_ps(_mm_dp_ps2(r12, r12, 0x7F), MINIMAL_CA_DISTANCE2)), 0)) {
          float e = ks_donor_acceptor(xyz, hcoords, nco_indices, ri, rj);
          if (e < HBOND_ENERGY_CUTOFF && !is_proline[ri])
            /* hbond from donor=ri to acceptor=rj */
            store_energies(hbonds, henergies, ri, rj, e);

          if (rj != ri + 1) {
            float e = ks_donor_acceptor(xyz, hcoords, nco_indices, rj, ri);
            if (e < HBOND_ENERGY_CUTOFF && !is_proline[rj])
              /* hbond from donor=rj to acceptor=ri */
              store_energies(hbonds, henergies, rj, ri, e);
          }
        }
      }
    }
    xyz += n_atoms*3; /* advance to the next frame */
    hbonds += n_residues*2;
    henergies += n_residues*2;
  }
  free(hcoords);
  free(skip);

  return 1;
}
