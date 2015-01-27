/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2012-2013 Stanford University and the Authors               */
/*                                                                       */
/* Authors: Robert McGibbon                                              */
/* Contributors:                                                         */
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
#include <stdio.h>
#include <math.h>
#include <pmmintrin.h>

#include "sasa.h"
#include "msvccompat.h"
#include "ssetools.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static void asa_frame(const float* frame, const int n_atoms, const float* atom_radii,
		      const float* sphere_points, const int n_sphere_points,
		      int* neighbor_indices, float* centered_sphere_points, float* areas)
{
  /*// Calculate the accessible surface area of each atom in a single snapshot
  //
  // Parameters
  // ----------
  // frame : 2d array, shape=[n_atoms, 3]
  //     The coordinates of the nuclei
  // n_atoms : int
  //     the major axis length of frame
  // atom_radii : 1d array, shape=[n_atoms]
  //     the van der waals radii of the atoms PLUS the probe radius
  // sphere_points : 2d array, shape=[n_sphere_points, 3]
  //     a bunch of uniformly distributed points on a sphere
  // n_sphere_points : int
  //    the number of sphere points

  // centered_sphere_points : WORK BUFFER 2d array, shape=[n_sphere_points, 3]
  //    empty memory that intermediate calculations can be stored in
  // neighbor_indices : WORK BUFFER 2d array, shape=[n_atoms]
  //    empty memory that intermediate calculations can be stored in
  // NOTE: the point of these work buffers is that if we want to call
  //    this function repreatedly, its more efficient not to keep re-mallocing
  //    these work buffers, but instead just reuse them.

  // areas : 1d array, shape=[n_atoms]
  //     the output buffer to place the results in -- the surface area of each
  //     atom
  */

  int i, j, k, k_prime;
  __m128 r, r_i, r_j, r_ij, atom_radius_i, atom_radius_j, radius_cutoff;
  __m128 radius_cutoff2, sp, r_jk, r2;
  int n_neighbor_indices, is_accessible, k_closest_neighbor;
  float constant = 4.0 * M_PI / n_sphere_points;

  for (i = 0; i < n_atoms; i++) {
    atom_radius_i = _mm_set1_ps(atom_radii[i]);
    r_i = load_float3(frame+i*3);

    /* Get all the atoms close to atom `i` */
    n_neighbor_indices = 0;
    for (j = 0; j < n_atoms; j++) {
      if (i == j) {
	continue;
      }

      r_j = load_float3(frame+j*3);
      r_ij = _mm_sub_ps(r_i, r_j);
      atom_radius_j = _mm_set1_ps(atom_radii[j]);

      /* Look for atoms `j` that are nearby atom `i` */
      radius_cutoff =  _mm_add_ps(atom_radius_i, atom_radius_j);
      radius_cutoff2 = _mm_mul_ps(radius_cutoff, radius_cutoff);
      r2 = _mm_dp_ps2(r_ij, r_ij, 0x7F);
      if (_mm_extract_epi16(CAST__M128I(_mm_cmplt_ps(r2, radius_cutoff2)), 0)) {
	    neighbor_indices[n_neighbor_indices]  = j;
	    n_neighbor_indices++;
      }
    if (_mm_extract_epi16(CAST__M128I(_mm_cmplt_ps(r2, _mm_set1_ps(1e-10))), 0)) {
	  printf("ERROR: THIS CODE IS KNOWN TO FAIL WHEN ATOMS ARE VIRTUALLY");
	  printf("ON TOP OF ONE ANOTHER. YOU SUPPLIED TWO ATOMS %f", _mm_cvtss_f32(r));
	  printf("APART. QUITTING NOW");
	  exit(1);
      }
    }

    /* Center the sphere points on atom i */
    for (j = 0; j < n_sphere_points; j++) {
      sp = _mm_add_ps(r_i, _mm_mul_ps(atom_radius_i, load_float3(sphere_points + 3*j)));
      store_float3(centered_sphere_points + 3*j, sp);
    }

    /* Check if each of these points is accessible */
    k_closest_neighbor = 0;
    for (j = 0; j < n_sphere_points; j++) {
      is_accessible = 1;
      r_j = load_float3(centered_sphere_points + 3*j);

      /* iterate through the sphere points by cycling through them */
      /* in a circle, starting with k_closest_neighbor and then wrapping */
      /* around */
      for (k = k_closest_neighbor; k < n_neighbor_indices + k_closest_neighbor; k++) {
	k_prime = k % n_neighbor_indices;
	r = _mm_set1_ps(atom_radii[neighbor_indices[k_prime]]);

	r_jk = _mm_sub_ps(r_j, load_float3(frame+3*neighbor_indices[k_prime]));
	if (_mm_extract_epi16(CAST__M128I(_mm_cmplt_ps(_mm_dp_ps2(r_jk, r_jk, 0xFF), _mm_mul_ps(r, r))), 0)) {
	  k_closest_neighbor = k;
	  is_accessible = 0;
	  break;
	}
      }

      if (is_accessible) {
	areas[i]++;
      }
    }

    areas[i] *= constant * (atom_radii[i])*(atom_radii[i]);
  }
}

static void generate_sphere_points(float* sphere_points, int n_points)
{
  /*
  // Compute the coordinates of points on a sphere using the
  // Golden Section Spiral algorithm.
  //
  // Parameters
  // ----------
  // sphere_points : array, shape=(n_points, 3)
  //     Empty array of length n_points*3 -- will be filled with the points
  //     as an array in C-order. i.e. sphere_points[3*i], sphere_points[3*i+1]
  //     and sphere_points[3*i+2] are the x,y,z coordinates of the ith point
  // n_pts : int
  //     Number of points to generate on the sphere
  //
  */
  int i;
  float y, r, phi;
  float inc = M_PI * (3.0 - sqrt(5));
  float offset = 2.0 / n_points;

  for (i = 0; i < n_points; i++) {
    y = i * offset - 1.0 + (offset / 2.0);
    r = sqrt(1.0 - y*y);
    phi = i * inc;

    sphere_points[3*i] = cos(phi) * r;
    sphere_points[3*i+1] = y;
    sphere_points[3*i+2] = sin(phi) * r;
  }
}


int sasa(const int n_frames, const int n_atoms, const float* xyzlist,
         const float* atom_radii, const int n_sphere_points,
         const int* atom_mapping, const int n_groups, float* out)
{
  /*
  // Calculate the accessible surface area of each atom in each frame of
  // a trajectory
  //
  // Parameters
  // ----------
  // xyzlist : 3d array, shape=[n_frames, n_atoms, 3]
  //     The coordinates of the nuclei
  // n_frames : int
  //     the number of frames in the trajectory
  // n_atoms : int
  //     the number of atoms in each frame
  // atom_radii : 1d array, shape=[n_atoms]
  //     the van der waals radii of the atoms PLUS the probe radius
  // n_sphere_points : int
  //     number of points to generate sampling the unit sphere. higher is
  //     better (more accurate) but more expensive
  // atom_mapping : 1d array, shape=[n_atoms]
  //     mapping from atoms onto groups, over which to accumulate the sasa.
  //     If `atom_mapping[i] = j`, that means that the ith atom is in group
  //     j. The groups must be contiguous integers starting from 0 to n_groups-1.
  // out : 2d array, shape=[n_frames, n_groups]
  //     the output buffer to place the results in. this array must be
  //     initialized with all zeros. out[i*n_groups + j] gives, in the `i`th frame
  //     of the trajectory, the total SASA of all of the atoms in group `j`.
  //
  */

  int i, j;

  /* work buffers that will be thread-local */
  int* wb1;
  float* wb2;
  float* outframe;
  float* outframebuffer;

  /* generate the sphere points */
  float* sphere_points = (float*) malloc(n_sphere_points*3*sizeof(float));
  generate_sphere_points(sphere_points, n_sphere_points);

#ifdef _OPENMP
  #pragma omp parallel private(wb1, wb2, outframebuffer, outframe)
  {
#endif

  /* malloc the work buffers for each thread */
  wb1 = (int*) malloc(n_atoms*sizeof(int));
  wb2 = (float*) malloc(3*n_sphere_points*sizeof(float));
  outframebuffer = (float*) calloc(n_atoms, sizeof(float));

#ifdef _OPENMP
  #pragma omp for
#endif
  for (i = 0; i < n_frames; i++) {
    asa_frame(xyzlist + i*n_atoms*3, n_atoms, atom_radii, sphere_points,
	      n_sphere_points, wb1, wb2, outframebuffer);
    outframe = out + (n_groups * i);
    for (j = 0; j < n_atoms; j++) {
        outframe[atom_mapping[j]] += outframebuffer[j];
    }
  }

  free(wb1);
  free(wb2);
  free(outframebuffer);
#ifdef _OPENMP
  } /* close omp parallel private */
#endif

  free(sphere_points);
  return 1;
}
