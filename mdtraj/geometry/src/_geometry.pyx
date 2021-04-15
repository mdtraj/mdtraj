# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import sys

##############################################################################
# Headers
##############################################################################

cdef extern from "geometry.h" nogil:
    void dist(const float* xyz, const int* pairs, float* distance_out,
              float* displacement_out, int n_frames,  int n_atoms,
              int n_pairs)
    void dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
                  float* distance_out, float* displacement_out,
                  int n_frames, int n_atoms, int n_pairs)
    void dist_t(const float* xyz, const int* pairs, const int* times,
                    float* distance_out, float* displacement_out, const int n_times, 
                    const int n_atoms, const int n_pairs)
    void dist_mic_t(const float* xyz, const int* pairs, const int* times,
                    const float* box_matrix, float* distance_out,
                    float* displacement_out, const int n_times, const int n_atoms,
                    const int n_pairs)
    void dist_mic_triclinic(const float* xyz, const int* pairs,
                            const float* box_matrix, float* distance_out,
                            float* displacement_out, int n_frames, int n_atoms,
                            int n_pairs)
    void dist_mic_triclinic_t(const float* xyz, const int* pairs, const int* times,
                              const float* box_matrix, float* distance_out,
                              float* displacement_out, int n_frames, int n_atoms,
                              int n_pairs)

    void angle(const float* xyz, const int* triplets, float* out,
               int n_frames, int n_atoms, int n_angles)

    void angle_mic(const float* xyz, const int* triplets, const float* box_matrix,
                   float* out, int n_frames, int n_atoms, int n_angles)

    void angle_mic_triclinic(const float* xyz, const int* triplets, const float* box_matrix,
                             float* out, int n_frames, int n_atoms, int n_angles)

    void dihedral(const float* xyz, const int* quartets, float* out,
                  int n_frames, int n_atoms,  int n_quartets)

    void dihedral_mic(const float* xyz, const int* quartets, float* out,
                      const float* box_matrix, int n_frames, int n_atoms,
                      int n_quartets)

    void dihedral_mic_triclinic(const float* xyz, const int* quartets, float* out,
                                const float* box_matrix, int n_frames, int n_atoms,
                                int n_quartets)

    void kabsch_sander(float* xyz, int* nco_indices, int* ca_indices,
                       const int* is_proline, int n_frames, int n_atoms,
                       int n_residues, int* hbonds, float* henergies)

    void dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
              const int* is_proline, const int* chains_ids, const int n_frames,
              const int n_atoms, const int n_residues, char* secondary)

    void find_closest_contact(const float* positions, const int* group1, const int* group2,
                              int n_group1, int n_group2, const float* box_vectors_pointer,
                              int* atom1, int* atom2, float* distance)

cdef extern from "sasa.h":
    void sasa(const int n_frames, const int n_atoms, const float* xyzlist,
              const float* atom_radii, const int n_sphere_points,
              const int* atom_mapping, const int n_groups, float* out) nogil


##############################################################################
# Wrappers
##############################################################################

def _dist(float[:, :, ::1] xyz,
          int[:, ::1] pairs,
          float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], &pairs[0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)


def _dist_displacement(float[:, :, ::1] xyz,
                       int[:, ::1] pairs,
                       float[:, :, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], <int*> &pairs[0,0], NULL, &out[0,0,0], n_frames, n_atoms, n_pairs)


def _dist_mic(float[:, :, ::1] xyz,
              int[:, ::1] pairs,
              float[:, :, ::1] box_matrix,
              float[:, ::1] out,
              orthogonal):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    if orthogonal:
        dist_mic(&xyz[0,0,0], &pairs[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)
    else:
        dist_mic_triclinic(&xyz[0,0,0], &pairs[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)


def _dist_t(float[:, :, ::1] xyz,
        int[:, ::1] pairs,
        int[:, ::1] times,
        float[:, ::1] out):
    cdef int n_times = times.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist_t(&xyz[0,0,0], &pairs[0,0], &times[0,0], &out[0,0], NULL, n_times, n_atoms, n_pairs)


def _dist_mic_t(float[:, :, ::1] xyz,
              int[:, ::1] pairs,
              int[:, ::1] times,
              float[:, :, ::1] box_matrix,
              float[:, ::1] out,
              orthogonal):
    cdef int n_times = times.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    if orthogonal:
        dist_mic_t(&xyz[0,0,0], &pairs[0,0], &times[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_times, n_atoms, n_pairs)
    else:
        dist_mic_triclinic_t(&xyz[0,0,0], &pairs[0,0], &times[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_times, n_atoms, n_pairs)


def _dist_mic_displacement(float[:, :, ::1] xyz,
                           int[:, ::1] pairs,
                           float[:, :, ::1] box_matrix,
                           float[:, :, ::1] out,
                           orthogonal):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    if orthogonal:
        dist_mic(&xyz[0,0,0], <int*> &pairs[0,0], &box_matrix[0,0,0], NULL, &out[0,0, 0], n_frames, n_atoms, n_pairs)
    else:
        dist_mic_triclinic(&xyz[0,0,0], <int*> &pairs[0,0], &box_matrix[0,0,0], NULL, &out[0,0, 0], n_frames, n_atoms, n_pairs)


def _angle(float[:, :, ::1] xyz,
           int[:, ::1] triplets,
           float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    angle(&xyz[0,0,0], &triplets[0,0], &out[0,0], n_frames, n_atoms, n_angles)


def _angle_mic(float[:, :, ::1] xyz,
               int[:, ::1] triplets,
               float[:, :, ::1] box_matrix,
               float[:, ::1] out,
               orthogonal):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    if orthogonal:
        angle_mic(&xyz[0,0,0], &triplets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_angles)
    else:
        angle_mic_triclinic(&xyz[0,0,0], &triplets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_angles)


def _dihedral(float[:, :, ::1] xyz,
              int[:, ::1] quartets,
              float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    dihedral(&xyz[0,0,0], <int*> &quartets[0,0], &out[0,0], n_frames, n_atoms, n_quartets)


def _dihedral_mic(float[:, :, ::1] xyz,
                  int[:, ::1] quartets,
                  float[:, :, ::1] box_matrix,
                  float[:, ::1] out,
                  orthogonal):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    if orthogonal:
        dihedral_mic(&xyz[0,0,0], <int*> &quartets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_quartets)
    else:
        dihedral_mic_triclinic(&xyz[0,0,0], <int*> &quartets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_quartets)


def _kabsch_sander(float[:, :, ::1] xyz,
                   int[:, ::1] nco_indices,
                   int[::1] ca_indices,
                   int[::1] is_proline,
                   int[:, :, ::1] hbonds,
                   float[:, :, ::1] henergies):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_residues = ca_indices.shape[0]
    kabsch_sander(&xyz[0,0,0], &nco_indices[0,0], &ca_indices[0],
                  &is_proline[0], n_frames, n_atoms, n_residues,
                  &hbonds[0,0,0], &henergies[0,0,0])


def _sasa(float[:, :, ::1] xyz,
          float[::1] atom_radii,
          int n_sphere_points,
          int[::1] atom_outmapping,
          float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    sasa(n_frames, n_atoms, &xyz[0,0,0], &atom_radii[0], n_sphere_points,
         &atom_outmapping[0], out.shape[1], &out[0,0])


def _dssp(float[:, :, ::1] xyz,
          int[:, ::1] nco_indices,
          int[::1] ca_indices,
          int[::1] is_proline,
          int[::1] chain_ids):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_residues = ca_indices.shape[0]
    cdef char[::1] secondary = bytearray(n_frames*n_residues)
    dssp(&xyz[0,0,0], &nco_indices[0,0], &ca_indices[0],
         &is_proline[0], &chain_ids[0], n_frames, n_atoms,
         n_residues, &secondary[0])

    PY2 = sys.version_info[0] == 2
    value = str(secondary.base) if PY2 else secondary.base.decode('ascii')
    return value


def _find_closest_contact(float[:, ::1] positions,
          int[::1] group1,
          int[::1] group2,
          box):
    cdef int atom1
    cdef int atom2
    cdef float distance
    cdef float[:, ::1] box_vectors
    cdef float* box_vectors_pointer
    if box is not None:
        box_vectors = np.asarray(box, order='c')
        box_vectors_pointer = &box_vectors[0,0]
    else:
        box_vectors_pointer = NULL
    find_closest_contact(&positions[0,0], &group1[0], &group2[0], len(group1), len(group2), box_vectors_pointer, &atom1, &atom2, &distance)
    return (atom1, atom2, distance)

include "image_molecules.pxi"
