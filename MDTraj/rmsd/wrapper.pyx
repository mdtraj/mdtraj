from cpython cimport bool
from mdtraj.utils.arrays import ensure_type
import cython
import numpy as np
cimport numpy as np
np.import_array()
from cython.parallel cimport prange

cdef extern float msd_axis_major(int nrealatoms, int npaddedatoms, int rowstride,
                                 float* aT, float* bT, float G_a, float G_b) nogil

cdef extern float msd_atom_major(int nrealatoms, int npaddedatoms,  float* a,
                                 float* b, float G_a, float G_b) nogil

cdef extern from "math.h":
    float sqrtf(float x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleRMSDs_axis_major(
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz1 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz2 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g1 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g2 not None,
int n_atoms,
int frame,
bool parallel=True):
    """getMultipleRMSDs_axis_major(xyz1, xyz2, g1, g2, n_atoms, frame, parallel=True)
    
    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in axis-major orders

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
    rmsds: np.ndarray, shape=(len(xyz2)),
        RMSDS between xyz1[frame] and all of xyz2
    """
    
    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms_padded = xyz1.shape[2]
    cdef float msd

    assert xyz1.ndim == 3 and xyz2.ndim == 3 and xyz1.shape[1] == 3 and xyz2.shape[1] == 3
    if not ((xyz1.shape[2] % 4 == 0) & (xyz2.shape[2] % 4 == 0)):
        raise ValueError("Input arrays must have third dimension of 4*n, "
                         "found %d and %d." % (xyz1.shape[2], xyz2.shape[2]))
    if frame >= xyz1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz1.shape[0]))

    cdef np.ndarray[dtype=np.float32_t, ndim=1] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd = msd_axis_major(n_atoms, n_atoms_padded, n_atoms_padded,
                                   &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd = msd_axis_major(n_atoms, n_atoms_padded, n_atoms_padded,
                                   &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)

    return distances


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleRMSDs_atom_major(
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz1 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz2 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g1 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g2 not None,
int n_atoms,
int frame,
bool parallel=True):
    """getMultipleRMSDs_atom_major(xyz1, xyz2, g1, g2, n_atoms, frame, parallel=True)
    
    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in atom-major orders

    Parameters
    ----------
    xyz1 : np.ndarray, shape=(n_frames, n_atoms_padded, 3), dtype=float32
        Coordinates of reference frame.
    xyz2 : np.ndarray, shape=(n_frames, n_atoms_padded, 3), dtype=float32
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
    rmsds: np.ndarray, shape=(len(xyz2)),
        RMSDS between xyz1[frame] and all of xyz2
    """
    
    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms_padded = xyz1.shape[1]
    cdef float msd

    assert xyz1.ndim == 3 and xyz2.ndim == 3 and xyz1.shape[2] == 3 and xyz2.shape[2] == 3
    if not ((xyz1.shape[1] % 4 == 0) & (xyz2.shape[1] % 4 == 0)):
        raise ValueError("Input arrays must have middle dimension of 4*n, "
                         "found %d and %d." % (xyz1.shape[1], xyz2.shape[1]))
    if frame >= xyz1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz1.shape[0]))

    cdef np.ndarray[dtype=np.float32_t, ndim=1] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd = msd_atom_major(n_atoms, n_atoms_padded,  &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd = msd_atom_major(n_atoms, n_atoms_padded, &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)

    return distances
