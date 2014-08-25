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
import cython
import numpy as np
cimport numpy as np
np.import_array()
assert sizeof(int) == sizeof(np.int32_t)
assert sizeof(float) == sizeof(np.float32_t)

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

cdef extern int sasa(int n_frames, int n_atoms, const float* xyzlist,
                     const float* atom_radii, int n_sphere_points,
                     float* array_of_areas) nogil
cdef extern from "hardware.h":
    int processorSupportsSSE41()


##############################################################################
# Wrappers
##############################################################################

def _processor_supports_sse41():
    """Does the current processor support SSE4.1 instructions?"""
    return processorSupportsSSE41()


@cython.boundscheck(False)
def _dist(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
          np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], <int*> &pairs[0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _dist_displacement(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
          np.ndarray[np.float32_t, ndim=3, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], <int*> &pairs[0,0], NULL, &out[0,0,0], n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _dist_mic(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
              np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
              np.ndarray[np.float32_t, ndim=3, mode='c'] box_matrix not None,
              np.ndarray[np.float32_t, ndim=2, mode='c'] out not None,
              int outflag=0):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist_mic(&xyz[0,0,0], <int*> &pairs[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _dist_mic_displacement(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
              np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
              np.ndarray[np.float32_t, ndim=3, mode='c'] box_matrix not None,
              np.ndarray[np.float32_t, ndim=3, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist_mic(&xyz[0,0,0], <int*> &pairs[0,0], &box_matrix[0,0,0], NULL, &out[0,0, 0], n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
def _angle(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
           np.ndarray[np.int32_t, ndim=2, mode='c'] triplets not None,
           np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    angle(&xyz[0,0,0], <int*> &triplets[0,0], &out[0,0], n_frames, n_atoms, n_angles)


@cython.boundscheck(False)
def _angle_mic(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
               np.ndarray[np.int32_t, ndim=2, mode='c'] triplets not None,
               np.ndarray[np.float32_t, ndim=3, mode='c'] box_matrix not None,
               np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    angle_mic(&xyz[0,0,0], <int*> &triplets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_angles)


@cython.boundscheck(False)
def _dihedral(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
              np.ndarray[np.int32_t, ndim=2, mode='c'] quartets not None,
              np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    dihedral(&xyz[0,0,0], <int*> &quartets[0,0], &out[0,0], n_frames, n_atoms, n_quartets)


@cython.boundscheck(False)
def _dihedral_mic(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
                  np.ndarray[np.int32_t, ndim=2, mode='c'] quartets not None,
                  np.ndarray[np.float32_t, ndim=3, mode='c'] box_matrix not None,
                  np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    dihedral_mic(&xyz[0,0,0], <int*> &quartets[0,0], &box_matrix[0,0,0], &out[0,0], n_frames, n_atoms, n_quartets)


@cython.boundscheck(False)
def _kabsch_sander(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
                   np.ndarray[np.int32_t, ndim=2, mode='c'] nco_indices not None,
                   np.ndarray[np.int32_t, ndim=1, mode='c'] ca_indices not None,
                   np.ndarray[np.int32_t, ndim=1, mode='c'] is_proline not None,
                   np.ndarray[np.int32_t, ndim=3, mode='c'] hbonds not None,
                   np.ndarray[np.float32_t, ndim=3, mode='c'] henergies not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_residues = ca_indices.shape[0]
    kabsch_sander(&xyz[0,0,0], <int*> &nco_indices[0,0], <int*> &ca_indices[0],
                  <int*> &is_proline[0], n_frames, n_atoms, n_residues,
                  <int*> &hbonds[0,0,0], &henergies[0,0,0])


@cython.boundscheck(False)
def _sasa(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.float32_t, ndim=1, mode='c'] atom_radii not None,
          int n_sphere_points,
          np.ndarray[np.float32_t, ndim=2, mode='c'] array_of_areas not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    sasa(n_frames, n_atoms, &xyz[0,0,0], &atom_radii[0], n_sphere_points, &array_of_areas[0,0])


@cython.boundscheck(False)
def _dssp(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.int32_t, ndim=2, mode='c'] nco_indices not None,
          np.ndarray[np.int32_t, ndim=1, mode='c'] ca_indices not None,
          np.ndarray[np.int32_t, ndim=1, mode='c'] is_proline not None,
          np.ndarray[np.int32_t, ndim=1, mode='c'] chain_ids not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_residues = ca_indices.shape[0]
    cdef char[:] secondary = bytearray(n_frames*n_residues)
    dssp(&xyz[0,0,0], <int*> &nco_indices[0,0], <int*> &ca_indices[0],
         <int*> &is_proline[0], <int*> &chain_ids[0], n_frames, n_atoms,
         n_residues, <char*> &secondary[0])

    PY2 = sys.version_info.major
    value = str(secondary.base) if PY2 else secondary.base.decode('ascii')
    return [value[i*n_residues:(i+1)*n_residues] for i in range(n_frames)]

