from cpython cimport bool
from mdtraj.utils.arrays import ensure_type
import cython
import numpy as np
cimport numpy as np
from cython.parallel cimport prange

cdef extern float ls_rmsd2_aligned_T_g (const int nrealatoms, const int npaddedatoms, const int rowstride,
                          const float* aT, const float* bT, const float G_a, const float G_b) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
def rmsd_one_to_all(
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz1 not None, 
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz2 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g1 not None, 
np.ndarray[np.float32_t, ndim=1, mode="c"] g2 not None, 
int n_atoms,
int frame,
bool parallel=True
):
    """Calculate the RMSD of several frames to a single frame.

    Parameters
    ----------
    xyz1 : np.ndarray, shape=(n_frames, 3, n_atoms_padded), dtype=float32
        Coordinates of reference frame.
    xyz2 : np.ndarray, shape=(n_frames, 3, n_atoms_padded), dtype=float32
        Coordinates of frames to iterate over
    g1 : np.ndarray, shape = (n_frames), dtype='float32'
        Pre-calculated G factors (traces) for each frame in xyz1
    g2 : np.ndarray, shape = (n_frames), dtype='float32'
        Pre-calculated G factors (traces) for each frame in xyz2
    n_atoms : int
        The number of atoms in the system.
    frame : int
        Index of the desired reference frame in xyz1.
    parallel : bool, default True
        If True, use openmp parallelization.

    Returns
    -------
    ndarray of distances of length len(xyz2)
    """
    
    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms_padded = xyz1.shape[2]
    cdef int true_stride = n_atoms_padded * 3
    
    ensure_type(val=xyz1, dtype='float32', ndim=3, name="xyz1", shape=(n_frames, 3, None))
    ensure_type(val=xyz2, dtype='float32', ndim=3, name="xyz2", shape=(n_frames, 3, None))
    if not ((xyz1.shape[2] % 4 == 0) & (xyz2.shape[2] % 4 == 0)):
        raise(ValueError("Input arrays must have third dimension of 4*n, found %d and %d." % (xyz1.shape[2], xyz2.shape[2])))
    
    cdef np.ndarray[np.float32_t, ndim=1] distances = np.zeros(n_frames, dtype='float32')

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            distances[i] = ls_rmsd2_aligned_T_g(n_atoms, n_atoms_padded, n_atoms_padded, &xyz1[frame,0,0], &xyz2[i,0,0], g1[frame], g2[i]) ** 0.5
    else:
        for i in range(n_frames, nogil=True):
            distances[i] = ls_rmsd2_aligned_T_g(n_atoms, n_atoms_padded, n_atoms_padded, &xyz1[frame,0,0], &xyz2[i,0,0], g1[frame], g2[i]) ** 0.5

    return distances
