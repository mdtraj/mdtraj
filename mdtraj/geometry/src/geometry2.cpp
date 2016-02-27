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
#include <float.h>

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
