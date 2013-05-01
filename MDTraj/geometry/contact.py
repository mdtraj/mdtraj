import numpy as np


def _compute_atom_distances_xyz(xyz, atom_pairs):
    """Compute the bond angles for each frame in xyz

    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float
        The cartesian coordinates
    atom_pairs : np.ndarray, shape[num_pairs, 2], dtype=int
        Each row gives the indices of two atoms.

    Returns
    -------
    distances : np.ndarray, shape=[n_frames, num_pairs], dtype=float

    Notes
    -----
    This is a reference implementation in python/numpy.
    """
    delta = np.diff(xyz[:, atom_pairs], axis=2)[:, :, 0]

    distances = (delta ** 2.).sum(-1) ** 0.5

    return distances


def compute_atom_distances(traj, atom_pairs):
    """Compute the bond angles for each frame in traj

    Parameters
    ----------
    tray : Trajectory
        Trajectory to compute distances in
    atom_pairs : np.ndarray, shape[num_pairs, 2], dtype=int
        Each row gives the indices of two atoms.

    Returns
    -------
    distances : np.ndarray, shape=[n_frames, num_pairs], dtype=float

    Notes
    -----
    This is a reference implementation in python/numpy.
    """
    return _compute_atom_distances_xyz(traj.xyz, atom_pairs)
