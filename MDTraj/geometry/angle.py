import numpy as np

def bond_angles_xyz(xyz, angle_indices):
    """Compute the bond angles for each frame in xyz
    
    This is a reference single threaded implementation in python/numpy
    
    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float
        The cartesian coordinates
    angle_indices : np.ndarray, shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle
        
    Returns
    -------
    angles : np.ndarray, shape=[n_frames, n_angles], dtype=float
    """
    
    n_frames = xyz.shape[0]
    angles = np.zeros((n_frames, len(angle_indices)))

    for i in xrange(n_frames):
        for j, (m, o, n) in enumerate(angle_indices):
            u_prime = xyz[i, m, :] - xyz[i, o, :]
            v_prime = xyz[i, n, :] - xyz[i, o, :]
            u_norm = np.linalg.norm(u_prime)
            v_norm = np.linalg.norm(v_prime)

            angles[i, j] = np.arccos(np.dot(u_prime, v_prime) /
                (u_norm * v_norm))

    return angles

def bond_angles(traj, angle_indices):
    """Compute the bond angles for each frame in traj
    
    This is a reference single threaded implementation in python/numpy
    
    Parameters
    ----------
    xyz : np.ndarray, shape=[n_frames, n_atoms, 3], dtype=float
        The cartesian coordinates
    angle_indices : np.ndarray, shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle
        
    Returns
    -------
    angles : np.ndarray, shape=[n_frames, n_angles], dtype=float
    """    
    return bond_angles_xyz(traj.xyz, angle_indices)
