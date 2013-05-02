import numpy as np
import IRMSD

def calculate_G(xyz):
    xyz = xyz.astype('float64')
    return (xyz ** 2.0).sum(-1).sum(-1)

class RMSDTrajectory():
    def __init__(self, traj):
        self.n_frames = traj.n_frames
        self.n_atoms = traj.n_atoms
        self.n_atoms_padded = 4 + traj.n_atoms - traj.n_atoms % 4
        
        traj.center_coordinates()

        self.xyz = np.zeros((traj.n_frames, 3, self.n_atoms_padded), dtype='float32')
        self.xyz[:,:,:self.n_atoms] = traj.xyz.transpose((0,2,1))
        self.G = (traj.xyz.astype('float64') ** 2.).sum(-1).sum(-1).astype('float32')

class RMSD():
    def __init__(self):
        pass


    def one_to_all(self, r_traj1, r_traj2, index1):
        """Calculate a vector of distances from one frame of the first trajectory
        to all of the frames in the second trajectory

        The distances calculated are from the `index1`th frame of `r_traj1`
        to the frames in `r_traj2`

        Parameters
        ----------
        r_traj : rmsd.TheoData
            Calculate RMSD to the index1 frame of this.
        r_traj_many : rmsd.TheoData
            Calculate RMSD to all of these frames.
        index1 : int
            index of reference frame in `r_traj`

        Returns
        -------
        Vector of distances of length len(r_traj_many)

        Notes
        -----
        If the omp_parallel optional argument is True, we use shared-memory
        parallelization in C to do this faster.
        """

        return IRMSD.rmsd_one_to_all(r_traj1.xyz, r_traj2.xyz, r_traj1.G, r_traj2.G, r_traj1.n_atoms, index1)
