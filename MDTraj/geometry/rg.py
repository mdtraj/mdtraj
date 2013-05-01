import numpy as np


def _compute_rg_xyz(xyz, masses=None):
    """Compute the Rg for every frame.

    Parameters
    ----------
    xyz : ndarray
        xyz coordinates
    masses : ndarray
        Transition matrix

    Returns
    -------
    Rg : ndarray
        Rg for every frame

    Notes
    -----
    If masses are none, assumes equal masses.
    """
    traj_length, num_atoms, num_dims = xyz.shape
    if not num_dims == 3:
        raise ValueError("What did you pass me?")
    if not xyz.dtype == np.float32:
        xyz = np.float32(xyz)
    if masses is None:
        masses = np.ones(num_atoms)

    weights = masses / masses.sum()

    mu = xyz.mean(1)
    centered = (xyz.transpose((1, 0, 2)) - mu).transpose((1, 0, 2))
    squared_dists = (centered ** 2).sum(2)
    Rg = (squared_dists * weights).sum(1) ** 0.5

    return Rg


def compute_rg(traj, masses=None):
    """Compute the Rg for every frame.

    Parameters
    ----------
    traj : Trajectory

    masses : ndarray, optional
        array of atom masses.

    Returns
    -------
    Rg : ndarray
        Rg for every frame

    Notes
    -----
    If masses are none, assumes equal masses.
    """
    return _compute_rg_xyz(traj.xyz, masses=masses)
