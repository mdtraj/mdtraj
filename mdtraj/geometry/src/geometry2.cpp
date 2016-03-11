/*=======================================================================*/
/* MDTraj: A Python Library for Loading, Saving, and Manipulating        */
/*         Molecular Dynamics Trajectories.                              */
/* Copyright 2012-2016 Stanford University and the Authors               */
/*                                                                       */
/* Authors: Peter Eastman                                                */
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

#include "geometry.h"
#include "vectorize_sse.h"
#include <cmath>
#include <cstddef>
#include <float.h>

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
//#include "anglekernels.h"
//#include "dihedralkernels.h"

#define COMPILE_WITH_PERIODIC_BOUNDARY_CONDITIONS
#include "distancekernels.h"
//#include "anglekernels.h"
//#include "dihedralkernels.h"

void dist_mic_triclinic(const float* xyz, const int* pairs, const float* box_matrix,
                        float* distance_out, float* displacement_out, const int n_frames,
                        const int n_atoms, const int n_pairs) {
    bool store_displacement = (displacement_out != NULL);
    bool store_distance = (distance_out != NULL);
    for (int i = 0; i < n_frames; i++) {
        // Load the periodic box vectors and make sure they're in reduced form.

        fvec4 box_vec1(box_matrix[0], box_matrix[3], box_matrix[6], 0);
        fvec4 box_vec2(box_matrix[1], box_matrix[4], box_matrix[7], 0);
        fvec4 box_vec3(box_matrix[2], box_matrix[5], box_matrix[8], 0);
        box_vec3 -= box_vec2*roundf(box_vec3[1]/box_vec2[1]);
        box_vec3 -= box_vec1*roundf(box_vec3[0]/box_vec1[0]);
        box_vec2 -= box_vec1*roundf(box_vec2[0]/box_vec1[0]);
        float recip_box_size[3] = {1.0f/box_vec1[0], 1.0f/box_vec2[1], 1.0f/box_vec3[2]};
        for (int j = 0; j < n_pairs; j++) {
            // Compute the displacement.

            fvec4 pos1(xyz + 3*pairs[2*j + 0]);
            fvec4 pos2(xyz + 3*pairs[2*j + 1]);
            fvec4 r12 = pos2-pos1;
            r12 -= box_vec3*round(r12[2]*recip_box_size[2]);
            r12 -= box_vec2*round(r12[1]*recip_box_size[1]);
            r12 -= box_vec1*round(r12[0]*recip_box_size[0]);

            // We need to consider 27 possible periodic copies.

            float min_dist2 = FLT_MAX;
            fvec4 min_r = r12;
            for (int x = -1; x < 2; x++) {
                fvec4 ra = r12 + box_vec1*x;
                for (int y = -1; y < 2; y++) {
                    fvec4 rb = ra + box_vec2*y;
                    for (int z = -1; z < 2; z++) {
                        fvec4 rc = rb + box_vec3*z;
                        float dist2 = dot3(rc, rc);
                        if (dist2 <= min_dist2) {
                            min_dist2 = dist2;
                            min_r = rc;
                        }
                    }
                }
            }

            // Store results.

            if (store_displacement) {
                float temp[4];
                min_r.store(temp);
                *displacement_out = temp[0];
                displacement_out++;
                *displacement_out = temp[1];
                displacement_out++;
                *displacement_out = temp[2];
                displacement_out++;
            }
            if (store_distance) {
                *distance_out = sqrtf(min_dist2);
                distance_out++;
            }
        }

        // Advance to the next frame.

        xyz += n_atoms*3;
        box_matrix += 9;
    }
}

/**
 * Identify the closest contact between two groups of atoms.
 */
void find_closest_contact(const float* positions, const int* group1, const int* group2, int n_group1, int n_group2,
                          const float* box_vectors_pointer, int* atom1, int* atom2, float* distance) {
    *atom1 = 0;
    *atom2 = 0;
    float distance2 = FLT_MAX;
    fvec4 box_vec1, box_vec2, box_vec3;
    float recip_box_size[3];
    if (box_vectors_pointer != NULL) {
        box_vec1 = fvec4(box_vectors_pointer);
        box_vec2 = fvec4(box_vectors_pointer+3);
        box_vec3 = fvec4(box_vectors_pointer+6);
        recip_box_size[0] = 1.0f/box_vectors_pointer[0];
        recip_box_size[1] = 1.0f/box_vectors_pointer[4];
        recip_box_size[2] = 1.0f/box_vectors_pointer[8];
    }
    for (int i = 0; i < n_group1; i++) {
        fvec4 pos1(positions + 3*group1[i]);
        for (int j = 0; j < n_group2; j++) {
            fvec4 pos2(positions + 3*group2[j]);
            fvec4 delta = pos1-pos2;
            if (box_vectors_pointer != NULL) {
                delta -= box_vec3*floorf(delta[2]*recip_box_size[2]+0.5f);
                delta -= box_vec2*floorf(delta[1]*recip_box_size[1]+0.5f);
                delta -= box_vec1*floorf(delta[0]*recip_box_size[0]+0.5f);
            }
            float r2 = dot3(delta, delta);
            if (r2 < distance2) {
                *atom1 = group1[i];
                *atom2 = group2[j];
                distance2 = r2;
            }
        }
    }
    *distance = sqrtf(distance2);
}
