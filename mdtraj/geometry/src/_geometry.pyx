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
import warnings
import cython
import numpy as np

##############################################################################
# Headers
##############################################################################

cdef extern from "geometry.h":
    int dist(const float* xyz, const int* pairs, float* distance_out,
                     float* displacement_out, int n_frames,  int n_atoms,
                     int n_pairs) nogil
    int dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
                         float* distance_out, float* displacement_out,
                         int n_frames, int n_atoms, int n_pairs) nogil

    int angle(const float* xyz, const int* triplets, float* out,
              int n_frames, int n_atoms, int n_angles) nogil

    int angle_mic(const float* xyz, const int* triplets, const float* box_matrix,
                  float* out, int n_frames, int n_atoms, int n_angles) nogil

    int dihedral(const float* xyz, const int* quartets, float* out,
                 int n_frames, int n_atoms,  int n_quartets) nogil

    int dihedral_mic(const float* xyz, const int* quartets, float* out,
                     const float* box_matrix, int n_frames, int n_atoms,
                     int n_quartets) nogil

    int kabsch_sander(float* xyz, int* nco_indices, int* ca_indices,
                      const int* is_proline, int n_frames, int n_atoms,
                      int n_residues, int* hbonds, float* henergies) nogil

    int dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
             const int* is_proline, const int* chains_ids, const int n_frames,
             const int n_atoms, const int n_residues, char* secondary) nogil

cdef extern from "sasa.h":
    int sasa(const int n_frames, const int n_atoms, const float* xyzlist,
             const float* atom_radii, const int n_sphere_points,
             const int* atom_mapping, const int n_groups, float* out) nogil

cdef extern from "hardware.h":
    int processorSupportsSSE41()

##############################################################################
# Wrappers
##############################################################################

def _processor_supports_sse41():
    """Does the current processor support SSE4.1 instructions?"""
    warnings.warn(('_processor_supports_sse41 is depicrated, and will be '
                  'removed in MDTraj 1.4'), DeprecationWarning)
    return processorSupportsSSE41()

@cython.boundscheck(False)
def _dist(float[:, :, ::1] xyz,
          int[:, ::1] pairs,
          float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], &pairs[0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _dist_displacement(float[:, :, ::1] xyz,
                       int[:, ::1] pairs,
                       float[:, :, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], <int*> &pairs[0,0], NULL, &out[0,0,0], n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _dist_mic(float[:, :, ::1] xyz,
              int[:, ::1] pairs,
              float[:, :, ::1] box_matrix,
              float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist_mic(&xyz[0,0,0], &pairs[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _dist_mic_displacement(float[:, :, ::1] xyz,
                           int[:, ::1] pairs,
                           float[:, :, ::1] box_matrix,
                           float[:, :, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist_mic(&xyz[0,0,0], <int*> &pairs[0,0], &box_matrix[0,0,0], NULL, &out[0,0, 0], n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _angle(float[:, :, ::1] xyz,
           int[:, ::1] triplets,
           float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    angle(&xyz[0,0,0], &triplets[0,0], &out[0,0], n_frames, n_atoms, n_angles)


@cython.boundscheck(False)
def _angle_mic(float[:, :, ::1] xyz,
               int[:, ::1] triplets,
               float[:, :, ::1] box_matrix,
               float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    angle_mic(&xyz[0,0,0], &triplets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_angles)


@cython.boundscheck(False)
def _dihedral(float[:, :, ::1] xyz,
              int[:, ::1] quartets,
              float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    dihedral(&xyz[0,0,0], <int*> &quartets[0,0], &out[0,0], n_frames, n_atoms, n_quartets)


@cython.boundscheck(False)
def _dihedral_mic(float[:, :, ::1] xyz,
                  int[:, ::1] quartets,
                  float[:, :, ::1] box_matrix,
                  float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    dihedral_mic(&xyz[0,0,0], <int*> &quartets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_quartets)


@cython.boundscheck(False)
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


@cython.boundscheck(False)
def _sasa(float[:, :, ::1] xyz,
          float[::1] atom_radii,
          int n_sphere_points,
          int[::1] atom_outmapping,
          float[:, ::1] out):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    sasa(n_frames, n_atoms, &xyz[0,0,0], &atom_radii[0], n_sphere_points,
         &atom_outmapping[0], out.shape[1], &out[0,0])


@cython.boundscheck(False)
def _dssp(float[:, :, ::1] xyz,
          int[:, ::1] nco_indices,
          int[::1] ca_indices,
          int[::1] is_proline,
          int[::1] chain_ids):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_residues = ca_indices.shape[0]
    cdef char[:] secondary = bytearray(n_frames*n_residues)
    dssp(&xyz[0,0,0], &nco_indices[0,0], &ca_indices[0],
         &is_proline[0], &chain_ids[0], n_frames, n_atoms,
         n_residues, &secondary[0])

    PY2 = sys.version_info[0] == 2
    value = str(secondary.base) if PY2 else secondary.base.decode('ascii')
    return value
