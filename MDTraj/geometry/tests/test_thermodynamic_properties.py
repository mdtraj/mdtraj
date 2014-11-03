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

from __future__ import print_function
import itertools
import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif, raises

tip3p_charges = np.array([-0.834, 0.417, 0.417])

def test_density():
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))
    rho = md.geometry.density(traj)
    
    eq(float(rho.mean()), 991., decimal=-1)

def test_dipole_moments():
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))
    
    charges = np.tile(tip3p_charges, traj.n_residues)

    moments0 = md.geometry.dipole_moments(traj, charges)

    # Now we screw up the molecule wholeness referencing all distances relative to atom zero.
    atom_indices = np.array([np.zeros(traj.n_atoms), np.arange(traj.n_atoms)], dtype='int32').T  # E.g. [[0, 0], [0, 1], [0, 2], [0, 3]]...
    xyz = md.compute_displacements(traj, atom_indices, periodic=True)  # Define coordinates relative to atom 0, PBC corrected.
    traj.xyz = xyz

    moments1 = md.geometry.dipole_moments(traj, charges)

    eq(moments0, moments1, decimal=4)
        

def test_static_dielectric():
    traj = md.load(get_fn("tip3p_300K_1ATM.xtc"), top=get_fn("tip3p_300K_1ATM.pdb"))
    
    charges = np.tile(tip3p_charges, traj.n_residues)

    temperature = 300

    epsilon0 = md.geometry.static_dielectric(traj, charges, temperature)

    atom_indices = np.array([np.zeros(traj.n_atoms), np.arange(traj.n_atoms)], dtype='int32').T  # E.g. [[0, 0], [0, 1], [0, 2], [0, 3]]...
    xyz = md.compute_displacements(traj, atom_indices, periodic=True)  # Define coordinates relative to atom 0, PBC corrected.
    traj.xyz = xyz
    epsilon1 = md.geometry.static_dielectric(traj, charges, temperature)

    eq(epsilon0, epsilon1, decimal=3)
    eq(float(epsilon0), 78.0, decimal=-1)
        
