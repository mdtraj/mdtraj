import numpy as np

from geometry cimport compute_distances as _compute_distances
from geometry cimport compute_distances_orthorhombic as _compute_distances_orthorhombic
from geometry cimport compute_distances_triclinic as _compute_distances_triclinic
from geometry cimport compute_angles as _compute_angles
from geometry cimport compute_angles_orthorhombic as _compute_angles_orthorhombic
from geometry cimport compute_dihedrals as _compute_dihedrals
from geometry cimport compute_dihedrals_orthorhombic as _compute_dihedrals_orthorhombic


def compute_distances(float[:, :, ::1] xyz, int[:, ::1] pairs):
     cdef int n_frames = xyz.shape[0]
     cdef int n_atoms = xyz.shape[1]
     cdef int n_pairs = pairs.shape[0]
     cdef float[:, ::1] out = np.zeros((n_frames, n_pairs), dtype=np.float32)

     _compute_distances(&xyz[0,0,0], &pairs[0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)
     return np.array(out, copy=False)


def compute_distances_orthorhombic(float[:, :, ::1] xyz, float[:, :, ::1] box_matrix, int[:, ::1] pairs):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    cdef float[:, ::1] out = np.zeros((n_frames, n_pairs), dtype=np.float32)

    _compute_distances_orthorhombic(&xyz[0,0,0], &pairs[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)
    return np.array(out, copy=False)


def compute_distances_triclinic(float[:, :, ::1] xyz, float[:, :, ::1] box_matrix, int[:, ::1] pairs):
    cdef int n_frames = xyz.shape[0]
    cdef int n_atoms = xyz.shape[1]
    cdef int n_pairs = pairs.shape[0]
    cdef float[:, ::1] out = np.zeros((n_frames, n_pairs), dtype=np.float32)

    _compute_distances_triclinic(&xyz[0,0,0], &pairs[0,0], &box_matrix[0,0,0], &out[0,0], NULL, n_frames, n_atoms, n_pairs)
    return np.array(out, copy=False)
