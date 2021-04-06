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
#include "vectorize.h"
#include <cstddef>
#include <math.h>
#include <float.h>
#include <vector>
#include "math_patch.h"

using std::vector;

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

#define COMPILE_WITH_TRICLINIC
#include "anglekernels.h"
#include "dihedralkernels.h"

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

            int offset1 = 3*pairs[2*j + 0];
            fvec4 pos1(xyz[offset1], xyz[offset1+1], xyz[offset1+2], 0);
            int offset2 = 3*pairs[2*j + 1];
            fvec4 pos2(xyz[offset2], xyz[offset2+1], xyz[offset2+2], 0);
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

void dist_mic_triclinic_t(const float* xyz, const int* pairs, const int* times,
                        const float* box_matrix, float* distance_out, 
                        float* displacement_out, const int n_frames,
                        const int n_atoms, const int n_pairs) {
    bool store_displacement = (displacement_out != NULL);
    bool store_distance = (distance_out != NULL);
    for (int i = 0; i < n_frames; i++) {
        // Load the periodic box vectors and make sure they're in reduced form.
        // Get box from frame of first index of time pair
        int box_offset = times[2*i + 0] * 9;
        box_matrix += box_offset;

        fvec4 box_vec1(box_matrix[0], box_matrix[3], box_matrix[6], 0);
        fvec4 box_vec2(box_matrix[1], box_matrix[4], box_matrix[7], 0);
        fvec4 box_vec3(box_matrix[2], box_matrix[5], box_matrix[8], 0);
        box_vec3 -= box_vec2*roundf(box_vec3[1]/box_vec2[1]);
        box_vec3 -= box_vec1*roundf(box_vec3[0]/box_vec1[0]);
        box_vec2 -= box_vec1*roundf(box_vec2[0]/box_vec1[0]);
        float recip_box_size[3] = {1.0f/box_vec1[0], 1.0f/box_vec2[1], 1.0f/box_vec3[2]};
        for (int j = 0; j < n_pairs; j++) {
            // Compute the displacement.

            int time_offset1 = 3*n_atoms*times[2*i + 0];
            int time_offset2 = 3*n_atoms*times[2*i + 1];
            int pair_offset1 = 3*pairs[2*j + 0];
            int pair_offset2 = 3*pairs[2*j + 1];
            int offset1 = time_offset1 + pair_offset1;
            int offset2 = time_offset2 + pair_offset2;
            fvec4 pos1(xyz[offset1], xyz[offset1+1], xyz[offset1+2], 0);
            fvec4 pos2(xyz[offset2], xyz[offset2+1], xyz[offset2+2], 0);
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
        // Reset box offset
        box_matrix -= box_offset;
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
        box_vec1 = fvec4(box_vectors_pointer[0], box_vectors_pointer[1], box_vectors_pointer[2], 0);
        box_vec2 = fvec4(box_vectors_pointer[3], box_vectors_pointer[4], box_vectors_pointer[5], 0);
        box_vec3 = fvec4(box_vectors_pointer[6], box_vectors_pointer[7], box_vectors_pointer[8], 0);
        recip_box_size[0] = 1.0f/box_vectors_pointer[0];
        recip_box_size[1] = 1.0f/box_vectors_pointer[4];
        recip_box_size[2] = 1.0f/box_vectors_pointer[8];
    }
    for (int i = 0; i < n_group1; i++) {
        int offset1 = 3*group1[i];
        fvec4 pos1(positions[offset1], positions[offset1+1], positions[offset1+2], 0);
        for (int j = 0; j < n_group2; j++) {
            int offset2 = 3*group2[j];
            fvec4 pos2(positions[offset2], positions[offset2+1], positions[offset2+2], 0);
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

/**
 * Conpute the Kabsch-Sander hydrogen bond energy between two residues
 * in a single conformation.
 */
static float ks_donor_acceptor(const float* xyz, const float* hcoords, const int* nco_indices, int donor, int acceptor) {
    fvec4 coupling(-2.7888f, -2.7888f, 2.7888f, 2.7888f); // 332 (kcal*A/mol) * 0.42 * 0.2 * (1nm / 10 A)
    fvec4 r_n(xyz[3*nco_indices[3*donor]], xyz[3*nco_indices[3*donor]+1], xyz[3*nco_indices[3*donor]+2], 0);
    fvec4 r_h(hcoords[4*donor], hcoords[4*donor+1], hcoords[4*donor+2], 0);
    fvec4 r_c(xyz[3*nco_indices[3*acceptor+1]], xyz[3*nco_indices[3*acceptor+1]+1], xyz[3*nco_indices[3*acceptor+1]+2], 0);
    fvec4 r_o(xyz[3*nco_indices[3*acceptor+2]], xyz[3*nco_indices[3*acceptor+2]+1], xyz[3*nco_indices[3*acceptor+2]+2], 0);
    fvec4 r_ho = r_h-r_o;
    fvec4 r_hc = r_h-r_c;
    fvec4 r_nc = r_n-r_c;
    fvec4 r_no = r_n-r_o;

    // Compute all four dot products (each of the squared distances) and pack them into a single fvec4.

    fvec4 d2_honchcno(dot3(r_ho, r_ho), dot3(r_nc, r_nc), dot3(r_hc, r_hc), dot3(r_no, r_no));
    fvec4 recip_sqrt = 1.0f/sqrt(d2_honchcno);
    float energy = dot4(coupling, recip_sqrt);
    return (energy < -9.9f ? -9.9f : energy);
}

/**
 * Assign hydrogen atom coordinates
 */
static void ks_assign_hydrogens(const float* xyz, const int* nco_indices, const int n_residues, float *hcoords, int* skip) {
    fvec4 r_n(xyz[3*nco_indices[0]], xyz[3*nco_indices[0]+1], xyz[3*nco_indices[0]+2], 0);
    r_n.store(hcoords);
    hcoords += 4;

    for (int ri = 1; ri < n_residues; ri++) {
        if (!skip[ri]) {
            int pc_index = nco_indices[3*(ri-1) + 1];
            int po_index = nco_indices[3*(ri-1) + 2];
            fvec4 pc(xyz[3*pc_index], xyz[3*pc_index+1], xyz[3*pc_index+2], 0);
            fvec4 po(xyz[3*po_index], xyz[3*po_index+1], xyz[3*po_index+2], 0);
            fvec4 r_co = pc-po;
            fvec4 r_n(xyz[3*nco_indices[3*ri]], xyz[3*nco_indices[3*ri]+1], xyz[3*nco_indices[3*ri]+2], 0);
            fvec4 norm_r_co = r_co/sqrt(dot3(r_co, r_co));
            fvec4 r_h = r_n+norm_r_co*0.1f;
            r_h.store(hcoords);
        }
        hcoords += 4;
    }
}

/**
 * Store a computed hbond energy and the appropriate residue indices
 * in the output arrays. This function is called twice by kabsch_sander,
 * so it seemed appropriate to factor it out.
 */
static void store_energies(int* hbonds, float* henergies, int donor, int acceptor, float e) {
    float existing_e0 = henergies[2*donor];
    float existing_e1 = henergies[2*donor+1];

    if (isnan(existing_e0) || e < existing_e0) {
        // Copy over any info in #0 hbond to #1
        hbonds[2*donor+1] = hbonds[donor*2];
        henergies[2*donor+1] = existing_e0;
        hbonds[2*donor] = acceptor;
        henergies[2*donor] = e;
    }
    else if (isnan(existing_e1) || e < henergies[2*donor+1]) {
        hbonds[2*donor+1] = acceptor;
        henergies[2*donor+1] = e;
    }
}

/**
 * Find all of backbone hydrogen bonds between residues in each frame of a trajectory.
 *
 * Parameters
 * ----------
 * xyz : array, shape=(n_frames, n_atoms, 3)
 *     The cartesian coordinates of all of the atoms in each frame.
 * nco_indices : array, shape=(n_residues, 3)
 *     The indices of the backbone N, C, and O atoms for each residue.
 * ca_indices : array, shape=(n_residues,)
 *     The index of the CA atom of each residue.
 * is_proline : array, shape=(n_residue,)
 *     If a particular residue does not contain a CA atom, or you want to skip
 *     the residue for another reason, this value should evaluate to True.
 *
 * Returns
 * -------
 * hbonds : array, shape=(n_frames, n_residues, 2)
 *     This is a little tricky, so bear with me. This array gives the indices
 *     of the residues that each backbone hbond *donor* is engaged in an hbond
 *     with. For instance, the equality `bonds[i, j, 0] == k` is interpreted as
 *     "in frame i, residue j is donating its first hydrogen bond from residue
 *     k". `bonds[i, j, 1] == k` means that residue j is donating its second
 *     hydrogen bond from residue k. A negative value indicates that no such
 *     hbond exists.
 * henergies : array, shape=(n_frames, n_residues, 2)
 *     The semantics of this array run parallel to the hbonds array, but
 *     instead of giving the identity of the interaction partner, it gives
 *     the energy of the hbond. Only hbonds with energy below -0.5 kcal/mol
 *     are recorded.
 */
void kabsch_sander(const float* xyz, const int* nco_indices, const int* ca_indices,
                   const int* is_proline, const int n_frames, const int n_atoms,
                   const int n_residues, int* hbonds, float* henergies) {
    float HBOND_ENERGY_CUTOFF = -0.5;
    float MINIMAL_CA_DISTANCE2 = 0.81f;
    vector<float> hcoords(n_residues*4);
    vector<int> skip(n_residues);
    for (int i = 0; i < n_residues; i++)
        if ((nco_indices[i*3] == -1) || (nco_indices[i*3+1] == -1) ||
            (nco_indices[i*3+2] == -1) || ca_indices[i] == -1)
            skip[i] = 1;

    for (int i = 0; i < n_frames; i++) {
        ks_assign_hydrogens(xyz, nco_indices, n_residues, &hcoords[0], &skip[0]);
        for (int ri = 0; ri < n_residues; ri++) {
            if (skip[ri])
                continue;
            fvec4 ri_ca(xyz[3*ca_indices[ri]], xyz[3*ca_indices[ri]+1], xyz[3*ca_indices[ri]+2], 0);
            for (int rj = ri + 1; rj < n_residues; rj++) {
                if (skip[rj])
                    continue;
                fvec4 rj_ca(xyz[3*ca_indices[rj]], xyz[3*ca_indices[rj]+1], xyz[3*ca_indices[rj]+2], 0);

                // Check the ca distance before proceding
                
                fvec4 r12 = ri_ca-rj_ca;
                if (dot3(r12, r12) < MINIMAL_CA_DISTANCE2) {
                    float e = ks_donor_acceptor(xyz, &hcoords[0], nco_indices, ri, rj);
                    if (e < HBOND_ENERGY_CUTOFF && !is_proline[ri])
                        store_energies(hbonds, henergies, ri, rj, e); // hbond from donor=ri to acceptor=rj
                    if (rj != ri + 1) {
                        float e = ks_donor_acceptor(xyz, &hcoords[0], nco_indices, rj, ri);
                        if (e < HBOND_ENERGY_CUTOFF && !is_proline[rj])
                            store_energies(hbonds, henergies, rj, ri, e); // hbond from donor=rj to acceptor=ri
                    }
                }
            }
        }
        xyz += n_atoms*3; // advance to the next frame
        hbonds += n_residues*2;
        henergies += n_residues*2;
    }
}
