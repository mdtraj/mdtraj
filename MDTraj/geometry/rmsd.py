import numpy as np
from mdtraj import IRMSD


def calculate_G(xyz):
    return (xyz.astype('float64') ** 2.0).sum(-1).sum(-1).astype('float32')

def reshape_irmsd(xyz):
    n_frames, n_atoms = xyz.shape[0:2]
    n_atoms_padded = 4 + n_atoms - n_atoms % 4
    xyz_irmsd = np.zeros((n_frames, 3, n_atoms_padded), dtype='float32')
    xyz_irmsd[:, :, :n_atoms] = xyz.transpose((0, 2, 1))
    return xyz_irmsd, n_atoms_padded

class RMSDTrajectory():
    def __init__(self, traj):
        self.n_frames = traj.n_frames
        self.n_atoms = traj.n_atoms

        traj.center_coordinates()
        self.G = calculate_G(traj.xyz)
        self.xyz, self.n_atoms_padded = reshape_irmsd(traj.xyz)

class RMSD():
    def __init__(self):
        pass

    def one_to_all(self, r_traj, r_traj_all, index1, parallel=True):
        """Calculate a vector of distances from one frame of the first trajectory
        to all of the frames in the second trajectory

        The distances calculated are from the `index1`th frame of `r_traj1`
        to the frames in `r_traj2`

        Parameters
        ----------
        r_traj : RMSDTrajectory
            Calculate RMSD to the index1 frame of this.
        r_traj_all : RMSDTrajectory
            Calculate RMSD to all of these frames.
        index1 : int
            index of reference frame in `r_traj`
        parallel : bool, default True
            if True, use OpenMP parallelization.

        Returns
        -------
        Vector of distances of length len(r_traj_all)
        """

        return IRMSD.rmsd_one_to_all(r_traj.xyz, r_traj_all.xyz, r_traj.G, r_traj_all.G, r_traj.n_atoms, index1, parallel)
