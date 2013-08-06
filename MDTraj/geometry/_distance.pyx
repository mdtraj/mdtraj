import numpy as np
cimport numpy as np

cdef extern int dist(float*, np.int32_t*, float*, float*, int, int, int)
cdef extern int dist_mic(float* xyz, np.int32_t* atom_pairs, float* box_matrix,
                         float* distance_out, float* distance_out,
                         int n_frames, int n_atoms, int n_pairs)

def distance(np.ndarray[dtype=np.float32_t, ndim=3] xyz not None,
             np.ndarray[dtype=np.int32_t, ndim=2]  pairs not None,
             np.ndarray[dtype=np.float32_t, ndim=3] box_vectors=None,
             np.ndarray[dtype=np.float32_t, ndim=2] out=None):
    """Compute the distances between atoms in atom pair for each frame in xyz.

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float32
        The cartesian coordinates of each atom in each simulation snapshot.
    pairs : np.ndarray, shape=[n_pairs, 2], dtype=int, default=None
        Each row gives the indices of two atoms whose distance we calculate.
    box_vectors : np.ndarray, shape=(3, 3), dtype=float32, default=None
        Periodic box vectors in each frame. If supplied, the distance will
        be computed using the minimum image convention.
    out : np.ndarray, shape=[n_frames, n_pairs], dtype=float32, optional
        You may presupply the output buffer
    
    Returns
    -------
    out : np.ndarray, shape=[n_frames, n_pairs], dtype=float32
         The distance, in each frame, between each pair of atoms.
    """
    if not xyz.flags.c_contiguous:
        raise ValueError('xyz must be c contiguous')

    if out is None:
        out = np.empty((xyz.shape[0], pairs.shape[0]), dtype=np.float32)
    else:
        if not out.flags.c_contiguous:
            raise ValueError('out must be c contiguous')

    if box_vectors is None:
        dist(&xyz[0,0,0], &pairs[0,0], &out[0,0], NULL,
             xyz.shape[0], xyz.shape[1], pairs.shape[0])
    else:
        dist_mic(&xyz[0,0,0], &pairs[0,0], &box_vectors[0, 0, 0], &out[0,0],
                 NULL, xyz.shape[0], xyz.shape[1], pairs.shape[0])

    return out


def displacement(np.ndarray[dtype=np.float32_t, ndim=3] xyz not None,
                 np.ndarray[dtype=np.int32_t, ndim=2]  pairs not None,
                 np.ndarray[dtype=np.float32_t, ndim=3] box_vectors=None,
                 np.ndarray[dtype=np.float32_t, ndim=3] out=None):
    """Compute the displacement vectors between pairs of atoms in each frame of
    xyz.

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float32
        The cartesian coordinates of each atom in each simulation snapshot.
    atom_pairs : np.ndarray, shape=[n_pairs, 2], dtype=int
        Each row gives the indices of two atoms whose displacement we calculate.
    box_vectors : np.ndarray, shape=(3, 3), dtype=float32, default=None
        Periodic box vectors in each frame. If supplied, the distance will
        be computed using the minimum image convention.
    out : np.ndarray, shape=[n_frames, n_pairs, 3], dtype=float32, optional
        You may presupply the output buffer

    Returns
    -------
    out : np.ndarray, shape=[n_frames, n_pairs, 3], dtype=float32
         The displacememt vector, in each frame, between each pair of atoms.
    """
    if not xyz.flags.c_contiguous:
        raise ValueError('xyz must be c contiguous')

    if out is None:
        out = np.empty((xyz.shape[0], pairs.shape[0], 3), dtype=np.float32)
    else:
        if not out.flags.c_contiguous:
            raise ValueError('out must be c contiguous')
        
    if box_vectors is None:
        dist(&xyz[0,0,0], &pairs[0,0], NULL, &out[0,0, 0],
             xyz.shape[0], xyz.shape[1], pairs.shape[0])
    else:
        if not box_vectors.flags.c_contiguous:
             raise ValueError('box_vectors must be c contiguous')
        dist_mic(&xyz[0,0,0], &pairs[0,0], &box_vectors[0, 0, 0], NULL,
                 &out[0,0,0], xyz.shape[0], xyz.shape[1], pairs.shape[0])

    return out
