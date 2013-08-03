import numpy as np
cimport numpy as np

cdef extern int dist(float*, np.int32_t*, float*, int, int, int)
cdef extern int dist_orthorhombic(float*, np.int32_t*, float*, float*)

def distance(np.ndarray[dtype=np.float32_t, ndim=3] xyz not None,
             np.ndarray[dtype=np.int32_t, ndim=2]  pairs not None,
             np.ndarray[dtype=np.float32_t, ndim=1] box_lengths=None,
             np.ndarray[dtype=np.float32_t, ndim=2] out=None):
    """Compute the distances between atoms in atom pair for each frame in xyz.

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float32
        The cartesian coordinates of each atom in each simulation snapshot.
    pairs : np.ndarray, shape=[n_pairs, 2], dtype=int, default=None
        Each row gives the indices of two atoms whose distance we calculate.
    box_lengths : np.ndarray, shape=(3,), dtype=float32, default=None
        If box lengths is supplied, and is not note, the distances will
        be computed with respect to the minimum image convention for an 
        orthorhombic box with the given box lengths.
    Returns
    -------
    distances : np.ndarray, shape=[n_frames, n_pairs], dtype=float32
         The distance, in each frame, between each pair of atoms.
    """
    if out is None:
        out = np.empty((xyz.shape[0], pairs.shape[0]), dtype=np.float32)

    if box_lengths is None:
        dist(&xyz[0,0,0], &pairs[0,0], &out[0,0],
             xyz.shape[0], xyz.shape[1], pairs.shape[0])
    else:
        dist_orthorhombic(&xyz[0,0,0], &pairs[0,0], &box_lengths[0], &out[0,0])

    return out


def displacement(np.ndarray[dtype=np.float32_t, ndim=3] xyz not None,
                 np.ndarray[dtype=np.int32_t, ndim=2]   atom_pairs not None,
                 np.ndarray[dtype=np.float32_t, ndim=1] box_lengths):
    """Compute the displacement vectors between pairs of atoms in each frame of
    xyz.

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float32
        The cartesian coordinates of each atom in each simulation snapshot.
    atom_pairs : np.ndarray, shape=[n_pairs, 2], dtype=int
        Each row gives the indices of two atoms whose displacement we calculate.
    box_lengths : np.ndarray, shape=(3,), dtype=float32, default=None
        If box lengths is supplied, and is not note, the distances will
        be computed with respect to the minimum image convention for an 
        orthorhombic box with the given box lengths.
    Returns
    -------
    displacement : np.ndarray, shape=[n_frames, n_pairs, 3], dtype=float32
         The displacememt vector, in each frame, between each pair of atoms.
    """
    pass
