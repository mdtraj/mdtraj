import numpy as np

from geometry cimport compute_distances as _compute_distances
from geometry cimport compute_angles as _compute_angles
from geometry cimport compute_dihedrals as _compute_dihedrals


def compute_distances(float[:, :, ::1] xyz, float[:, :, ::1] box_matrix, ssize_t[:, ::1] pairs):
    cdef ssize_t n_frames = xyz.shape[0]
    cdef ssize_t n_atoms = xyz.shape[1]
    cdef ssize_t n_pairs = pairs.shape[0]
    cdef float[:, ::1] out = np.zeros((n_frames, n_pairs), dtype=np.float32)
    cdef float* box_pointer = &box_matrix[0,0,0] if box_matrix is not None else NULL


    _compute_distances(&xyz[0,0,0], &pairs[0,0], box_pointer, &out[0,0], NULL, n_frames, n_atoms, n_pairs)
    return np.array(out, copy=False)


def compute_angles(float[:, :, ::1] xyz, float[:, :, ::1] box_matrix, ssize_t[:, ::1] triplets):
    cdef ssize_t n_frames = xyz.shape[0]
    cdef ssize_t n_atoms = xyz.shape[1]
    cdef ssize_t n_triplets = triplets.shape[0]
    cdef float[:, ::1] out = np.zeros((n_frames, n_triplets), dtype=np.float32)
    cdef float* box_pointer = &box_matrix[0,0,0] if box_matrix is not None else NULL

    _compute_angles(&xyz[0,0,0], &triplets[0,0], box_pointer, &out[0,0], n_frames, n_atoms, n_triplets)
    return np.array(out, copy=False)


def compute_dihedrals(float[:, :, ::1] xyz, float[:, :, ::1] box_matrix, ssize_t[:, ::1] quartets):
    cdef ssize_t n_frames = xyz.shape[0]
    cdef ssize_t n_atoms = xyz.shape[1]
    cdef ssize_t n_quartets = quartets.shape[0]
    cdef float[:, ::1] out = np.zeros((n_frames, n_quartets), dtype=np.float32)
    cdef float* box_pointer = &box_matrix[0,0,0] if box_matrix is not None else NULL


    _compute_dihedrals(&xyz[0,0,0], &quartets[0,0], box_pointer, &out[0,0], n_frames, n_atoms, n_quartets)
    return np.array(out, copy=False)
