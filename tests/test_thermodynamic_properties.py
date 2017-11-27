##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import mdtraj as md
import numpy as np
import pandas as pd
import pytest
from mdtraj.testing import eq

temperature = 300.
tip3p_charges = np.array([-0.834, 0.417, 0.417])

"""
Reference value taken from Gromacs 5.0.4 using the following commands.
Note that gromacs accumulates property averages over every timestep, so precise
comparisons require writing XTC output very frequently.


mkdir gromacs
cd gromacs
ln -s ../MDTraj/testing/reference/ ./reference

pdb2gmx -ff amber99sb -water tip3p -f ./reference/tip3p_300K_1ATM.pdb  -ignh
grompp -f ./reference/md.mdp -c ./conf.gro -p ./topol.top

# Do a quick equilibration run
export OMP_NUM_THREADS=1
mdrun -nsteps 500000

# Do a short production run
grompp -f ./reference/md.mdp -c ./confout.gro -p ./topol.top
mdrun -nsteps 200000
cp traj_comp.xtc reference/tip3p_300K_1ATM.xtc

# Compute tab separated dataset for energies
g_energy -xvg None
15
cp energy.xvg reference/tip3p_300K_1ATM.csv
(save enthalpy to csv style output)

# Compute static dielectric constant
g_dipoles

Epsilon = 87.1818

# Compute fluctuation properties
g_energy -fluct_props

Energy                      Average   Err.Est.       RMSD  Tot-Drift
-------------------------------------------------------------------------------
Kinetic En.                 1916.68          3    66.4762    9.59656  (kJ/mol)
Total Energy               -8286.19        9.9    120.584    28.1681  (kJ/mol)
Temperature                 298.411       0.46    10.3498     1.4941  (K)
Pressure                   -10.4124         12    606.328   -47.4295  (bar)
Volume                      7.92813     0.0086  0.0819152  0.0164395  (nm^3)
Density                     973.651        1.1    10.0623   -2.02709  (kg/m^3)
Enthalpy                   -8285.72        9.9    120.586    28.1691  (kJ/mol)


Volume                                   = 0.00477443 m^3/mol
Enthalpy                                 =   -8285.72 kJ/mol
Coefficient of Thermal Expansion Alpha_P = 0.000895685 (1/K)
Isothermal Compressibility Kappa         = 2.05427e-10 (J/m^3)
Adiabatic bulk modulus                   = 4.86791e+09 (m^3/J)
Heat capacity at constant pressure Cp    =    19639.5 J/mol K
Cp-Cv                                    =     5563.55 J/mol K

"""


def test_volume(get_fn):
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))
    v = traj.unitcell_volumes.mean()
    reference = 7.92813  # From gromacs, see above comment

    assert abs((v - reference) / reference) < 1E-3, "Volume tolerance not met!"


def test_density(get_fn):
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))
    rho = md.geometry.density(traj).mean()
    reference = 973.651  # From gromacs, see above comment

    assert abs((rho - reference) / reference) < 1E-3, "Density tolerance not met!"


def test_dipole_moments(get_fn):
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))

    charges = np.tile(tip3p_charges, traj.n_residues)

    moments0 = md.geometry.dipole_moments(traj, charges)

    # Now we screw up the molecule wholeness referencing all distances relative to atom zero.
    # E.g. [[0, 0], [0, 1], [0, 2], [0, 3]]...
    atom_indices = np.array([np.zeros(traj.n_atoms), np.arange(traj.n_atoms)], dtype='int32').T
    # Define coordinates relative to atom 0, PBC corrected.
    xyz = md.compute_displacements(traj, atom_indices, periodic=True)
    traj.xyz = xyz

    moments1 = md.geometry.dipole_moments(traj, charges)

    eq(moments0, moments1, decimal=4)


def test_static_dielectric(get_fn):
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))

    charges = np.tile(tip3p_charges, traj.n_residues)

    epsilon0 = md.geometry.static_dielectric(traj, charges, temperature)
    # E.g. [[0, 0], [0, 1], [0, 2], [0, 3]]...
    atom_indices = np.array([np.zeros(traj.n_atoms), np.arange(traj.n_atoms)], dtype='int32').T
    # Define coordinates relative to atom 0, PBC corrected.
    xyz = md.compute_displacements(traj, atom_indices, periodic=True)
    traj.xyz = xyz
    epsilon1 = md.geometry.static_dielectric(traj, charges, temperature)

    eq(epsilon0, epsilon1, decimal=3)

    reference = 87.1818  # From gromacs, see above comment

    assert abs((epsilon1 - reference) / reference) < 1E-3, "Dielectric tolerance not met!"


def test_kappa(get_fn):
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))
    kappa = md.geometry.isothermal_compressability_kappa_T(traj, temperature)
    reference = 2.05427E-10 * 1E5  # m^3 / J to 1 / bar.  Data from gromacs.  See above comment

    # 20% tolerance.
    assert abs((kappa - reference) / reference) < 2E-1, "Compressability tolerance not met!"


@pytest.mark.skip("Thermal expansion test doesn't work.")  # Not working
def test_alpha(get_fn):
    # Had some issues finding a useful unit test, so thermal_expansion_alpha_P() is currently disabled.
    # Feel free to file a pull request with a working unit test :)
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))

    data = pd.read_table(get_fn('tip3p_300K_1ATM.tab'), names=["timestep", "energy"], sep=r"\s*")
    energy = data.energy.values
    alpha = md.geometry.thermal_expansion_alpha_P(traj, temperature, energy)

    reference = 0.000895685  # From gromacs, see notes above.md
    eq(alpha, reference, decimal=3)
