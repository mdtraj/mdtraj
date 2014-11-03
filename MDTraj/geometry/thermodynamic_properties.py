import numpy as np
import mdtraj as md
import mdtraj.utils.unit.unit_definitions as u
from mdtraj.utils.unit import in_units_of

# Units taken from http://en.wikipedia.org/wiki/Boltzmann_constant on Nov. 2.
kB = 1.3806488E-23 * u.joule / u.kelvin
epsilon0 = 8.854187817E-12 * u.farad / u.meter

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
    """Calculate the static dielectric constant from a trajectory.

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
    
    dipole_variance = (subtracted * subtracted).sum(-1).mean(0) * (u.elementary_charge * u.nanometers) ** 2.  # <M*M> - <M>*<M>
    
    volume = np.mean(map(np.linalg.det, traj.unitcell_vectors)) * u.nanometers ** 3.  # Average box volume of trajectory
    
    static_dielectric = 1.0 + dipole_variance / (3 * kB * temperature * volume * epsilon0)  # Eq. 7 of Derivation of an improved simple point charge model for liquid water: SPC/A and SPC/L 
    
    return static_dielectric


def heat_capacity_Cp():
    raise(NotImplementedError("This has not been implemented yet!"))

def isothermal_compressability_kappa_T():
    raise(NotImplementedError("This has not been implemented yet!"))

def thermal_expansion_alpha_P():
    raise(NotImplementedError("This has not been implemented yet!"))


def density(traj):
    """Calculate the instantaneous mass density of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mdtraj trajectory.

    Returns
    -------
    density_trace : np.array, shape=(n_frames), dtype=float
        A timeseries of the mass density of each frame.
    """
    mass = sum(atom.element.mass for atom in traj.top.atoms)
    volume_trace = np.array(map(np.linalg.det, traj.unitcell_vectors))
    densities = mass / volume_trace

    #conversion = 1.660538921E-23
    conversion *= in_units_of(1., "nanometer ** -3", "meter ** -3")
    conversion = in_units_of(1., "dalton * nanometer ** -3", "kilogram / item * meter ** -3")
    
    x = 1.0 * u.dalton
    rho = x / (600 * u.nanometer ** 3)

    rho.value_in_unit(u.gram / u.item / u.milliliter)


    return densities * conversion
