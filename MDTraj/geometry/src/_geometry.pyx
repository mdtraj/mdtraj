import cython
import numpy as np
cimport numpy as np
np.import_array()
assert sizeof(int) == sizeof(np.int32_t)
assert sizeof(float) == sizeof(np.float32_t)

cdef extern int dist(float* xyz, int* pairs, float* distance_out,
                     float* displacement_out, int n_frames,  int n_atoms,
                     int n_pairs) nogil

cdef extern int dist_mic(float* xyz, int* pairs, float* box_matrix,
                         float* distance_out, float* displacement_out,
                         int n_frames, int n_atoms, int n_pairs) nogil
                         
cdef extern int angle(float* xyz, int* triplets, float* out,
                      int n_frames, int n_atoms, int n_angles) nogil

cdef extern int dihedral(float* xyz, int* quartets, float* out,
                           int n_frames, int n_atoms,  int n_quartets) nogil
         
cdef extern int kabsch_sander(float* xyz, int* nco_indices, int* ca_indices,
                              int n_frames, int n_atoms, int n_residues,
                              int* hbonds, float* henergies) nogil

cdef extern int sasa(int n_frames, int n_atoms, float* xyzlist,
	                 float* atom_radii, int n_sphere_points, float* array_of_areas) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
def _dist(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
          np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], <int*> &pairs[0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)
    
@cython.boundscheck(False)
@cython.wraparound(False)
def _dist_displacement(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
          np.ndarray[np.float32_t, ndim=3, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist(&xyz[0,0,0], <int*> &pairs[0,0], NULL, &out[0,0,0], n_frames, n_atoms, n_pairs)



@cython.boundscheck(False)
@cython.wraparound(False)
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
@cython.wraparound(False)
def _dist_mic_displacement(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
              np.ndarray[np.int32_t, ndim=2, mode='c'] pairs not None,
              np.ndarray[np.float32_t, ndim=3, mode='c'] box_matrix not None,
              np.ndarray[np.float32_t, ndim=3, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    dist_mic(&xyz[0,0,0], <int*> &pairs[0,0], &box_matrix[0,0,0], NULL, &out[0,0, 0], n_frames, n_atoms, n_pairs)


@cython.boundscheck(False)
@cython.wraparound(False)
def _angle(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
           np.ndarray[np.int32_t, ndim=2, mode='c'] triplets not None,
           np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_angles = triplets.shape[0]
    print("CALLING ANGLE FROM CYTHON")
    angle(&xyz[0,0,0], <int*> &triplets[0,0], &out[0,0], n_frames, n_atoms, n_angles)
    print("RETURNED FROM C CALL")



@cython.boundscheck(False)
@cython.wraparound(False)
def _dihedral(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
              np.ndarray[np.int32_t, ndim=2, mode='c'] quartets not None,
              np.ndarray[np.float32_t, ndim=2, mode='c'] out not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_quartets = quartets.shape[0]
    dihedral(&xyz[0,0,0], <int*> &quartets[0,0], &out[0,0], n_frames, n_atoms, n_quartets)


@cython.boundscheck(False)
@cython.wraparound(False)
def _kabsch_sander(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
                   np.ndarray[np.int32_t, ndim=2, mode='c'] nco_indices not None,
                   np.ndarray[np.int32_t, ndim=2, mode='c'] ca_indices not None,
                   np.ndarray[np.int32_t, ndim=3, mode='c'] hbonds not None,
                   np.ndarray[np.float32_t, ndim=3, mode='c'] henergies not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_residues = ca_indices.shape[0]
    kabsch_sander(&xyz[0,0,0], <int*> &nco_indices[0,0], <int*> &ca_indices[0,0],
                  n_frames, n_atoms, n_residues, <int*> &hbonds[0,0,0], &henergies[0,0,0])


@cython.boundscheck(False)
@cython.wraparound(False)
def _sasa(np.ndarray[np.float32_t, ndim=3, mode='c'] xyz not None,
          np.ndarray[np.float32_t, ndim=1, mode='c'] atom_radii not None,
          int n_sphere_points,
          np.ndarray[np.float32_t, ndim=2, mode='c'] array_of_areas not None):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    sasa(n_frames, n_atoms, &xyz[0,0,0], &atom_radii[0], n_sphere_points, &array_of_areas[0,0])          
