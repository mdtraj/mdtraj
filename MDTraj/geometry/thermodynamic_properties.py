import numpy as np
import simtk.unit as u
import mdtraj as md

kB = u.BOLTZMANN_CONSTANT_kB

mu0 = 4 * np.pi * 1E-7 * u.henry / u.meter
epsilon0 = 1. / (mu0 * u.constants.SPEED_OF_LIGHT_C ** 2.)

def dipole_moments(traj, charges):
    """Calculate the dipole moments of each frame in a trajectory.


    Parameters
    ----------
    traj : Trajectory
        An mdtraj trajectory.
    charges : np.ndarray, shape=(n_atoms), dtype=float
       Charges of each atom in the topology, expressed in units of the
       elementary charge constant.  

    Returns
    -------
    moments : np.ndarray, shape=[n_frames, 3], dtype=float
        Dipole moments of trajectory, units of nm * elementary charge.
    """
    
    atom_indices = np.array([np.zeros(traj.n_atoms), np.arange(traj.n_atoms)]).T  # E.g. [[0, 0], [0, 1], [0, 2], [0, 3]]...
    
    xyz = md.compute_displacements(traj, atom_indices, periodic=True)  # Define coordinates relative to atom 0, PBC corrected.
    
    moments = xyz.transpose(0, 2, 1).dot(charges)  # Sum over atoms

    return moments

def static_dielectric(traj, charges, temperature):
    """Calculate the dipole moments of a trajectory.


    Parameters
    ----------
    traj : Trajectory
        An mdtraj trajectory.
    charges : np.ndarray, shape=(n_atoms), dtype=float
       Charges of each atom in the topology, expressed in units of the
       elementary charge constant.  
    temperature : temperature, simtk.unit.Quantity [Temperature]
        The temperature of interest.

    Returns
    -------
    static_dielectric : float, 
        The (unitless) relative static dielectric constant.
    """    
    moments = dipole_moments(traj, charges)
    
    mu = moments.mean(0)  # Mean over frames
    
    subtracted = moments - mu
    
    dipole_variance = (subtracted * subtracted).sum(-1).mean(0) * (u.elementary_charge * u.nanometers) ** 2.  # <M^2> - <M>^2
    
    volume = np.mean(map(np.linalg.det, traj.unitcell_vectors)) * u.nanometers ** 3.  # Average box volume of trajectory
    
    static_dielectric = 1.0 + dipole_variance * (4 * np.pi) / (3 * kB * temperature * volume * epsilon0)  # Eq. 2 of Fennell, 2012. J Phys Chem. B, except added epsilon0 in denominator.
    
    return static_dielectric


def heat_capacity_Cp():
    raise(NotImplementedError("This has not been implemented yet!"))

def isothermal_compressability_kappa_T():
    raise(NotImplementedError("This has not been implemented yet!"))

def thermal_expansion_alpha_P():
    raise(NotImplementedError("This has not been implemented yet!"))
