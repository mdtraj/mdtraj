"""
Notes
-----
The functions in this file use the MDTraj version of openmm.unit *internally*.
However, all inputs and outputs are done with implicit units.  This is
to avoid incompatibilities between versions of openmm.unit and MDTraj.utils.unit.

"""

import numpy as np

import mdtraj as md
import mdtraj.utils.unit.unit_definitions as u

# Units taken from http://en.wikipedia.org/wiki/Boltzmann_constant on Nov. 2.
kB = 1.3806488e-23 * u.joule / u.kelvin
epsilon0 = 8.854187817e-12 * u.farad / u.meter
gas_constant = 8.3144621 * u.joule / u.kelvin / u.mole


def dipole_moments(traj, charges):
    """Calculate the dipole moments of each frame in a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mdtraj trajectory.
    charges : np.ndarray, shape=(n_atoms), dtype=float
       Charges of each atom in the trajectory, expressed in units of the
       elementary charge constant.

    Returns
    -------
    moments : np.ndarray, shape=(n_frames, 3), dtype=float
        Dipole moments of trajectory, units of nm * elementary charge.

    Notes
    -----
    This code works by first calculating displacements relative to the
    first atom in each residue (e.g. local frame).  Then, this is added
    to the PBC-corrected displacement between the first atom in the two
    molecules.  This total displacement is then used as to calculate the
    box dipole moment.
    """
    local_indices = np.array(
        [(a.index, a.residue.atom(0).index) for a in traj.top.atoms],
        dtype="int32",
    )
    local_displacements = md.compute_displacements(traj, local_indices, periodic=True)

    molecule_indices = np.array(
        [(a.residue.atom(0).index, 0) for a in traj.top.atoms],
        dtype="int32",
    )
    molecule_displacements = md.compute_displacements(
        traj,
        molecule_indices,
        periodic=True,
    )

    xyz = local_displacements + molecule_displacements

    moments = xyz.transpose(0, 2, 1).dot(charges)

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
    temperature : temperature, float
        The temperature of interest, in units of kelvin.

    Returns
    -------
    static_dielectric : float,
        The (unitless) relative static dielectric constant.

    Notes
    -----
    See eqn. (2) in 10.1021/jp3002383 or eqn. (7) in 10.1063/1.1476316
    or https://github.com/gromacs/gromacs/blob/master/src/gromacs/gmxana/gmx_current.c#L622
    """
    temperature = temperature * u.kelvin

    moments = dipole_moments(traj, charges)

    mu = moments.mean(0)  # Mean over frames

    subtracted = moments - mu

    dipole_variance = (subtracted * subtracted).sum(-1).mean(0) * (
        u.elementary_charge * u.nanometers
    ) ** 2.0  # <M*M> - <M>*<M> = <(M - <M>) * (M - <M>)>

    volume = traj.unitcell_volumes.mean() * u.nanometers**3.0  # Average box volume of trajectory

    static_dielectric = 1.0 + dipole_variance / (
        3 * kB * temperature * volume * epsilon0
    )  # Eq. 7 of Derivation of an improved simple point charge model for liquid water: SPC/A and SPC/L
    # Also https://github.com/gromacs/gromacs/blob/master/src/gromacs/gmxana/gmx_current.c#L622

    return static_dielectric


def heat_capacity_Cp():
    raise (NotImplementedError("Has not been implemented."))


def isothermal_compressability_kappa_T(traj, temperature):
    """Calculate the isothermal compressability.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        An mdtraj trajectory.
    temperature : float
        The temperature of your trajectory, in units of kelvin.

    Returns
    -------
    kappa : float
        The isothermal compressability, in units of bar^-1.

    Notes
    -----
    Equation (4) in Fennell, Dill.  J. Phys. Chem. B, 2012.
    """
    temperature = temperature * u.kelvin

    mu = traj.unitcell_volumes.mean()

    kappa = np.cov(traj.unitcell_volumes) / mu
    kappa = kappa * u.nanometers**3

    kappa /= kB * temperature

    return kappa * u.bar


def thermal_expansion_alpha_P(traj, temperature, energies):
    """Calculate the thermal expansion coefficient.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        An mdtraj trajectory.
    temperature : float
        The temperature of your trajectory, in units of kelvin.
    energies : np.ndarray, dtype=float, shape=(n_frames)
        An array containing the potentail energies of each trajectory
        frame, in units of kJ / mol.

    Returns
    -------
    alpha : float
        The thermal expanssion coefficient, units of inverse Kelvin

    Notes
    -----
    Equation (5) in Fennell, Dill.  J. Phys. Chem. B, 2012.
    THIS FUNCTION IS NOT CURRENTLY IMPLEMENTED!
    """
    raise (NotImplementedError("Disabled due to lack of available unit test."))
    # Had some issues finding a useful unit test, so disabled this code for now.
    # Feel free to file a pull request with a working unit test :)
    temperature = temperature * u.kelvin

    mean_volume = traj.unitcell_volumes.mean()

    alpha = np.cov(traj.unitcell_volumes, energies)[0, 1]  # <HV> - <H><V> = cov(H, V)
    alpha /= mean_volume
    alpha *= u.kilojoules_per_mole

    alpha /= gas_constant * temperature**2

    return alpha * u.kelvin


def density(traj, masses=None):
    """Calculate the mass density of each frame in a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mdtraj trajectory.
    masses : np.ndarray, optional, default=None
        If not None, use these masses when calculating the density.  If
        None, then use the standard elemental masses associated with
        traj.topology.atoms.

    Returns
    -------
    density_trace : np.array, shape=(n_frames), dtype=float
        The mass density of each frame.

    """
    if masses is None:
        mass = sum(atom.element.mass for atom in traj.top.atoms)
    else:
        mass = sum(masses)

    volume_trace = traj.unitcell_volumes
    densities = mass / volume_trace

    # conversion = in_units_of(1., "dalton * nanometer ** -3", "kilogram / item * meter ** -3")
    # The item thing is really weird, but taken from OpenMM StateDataReporter's density calculation
    conversion = 1.6605387823355087  # The units stuff is busted on py3k, so using hard-coded for now.

    densities = densities * conversion
    return densities
