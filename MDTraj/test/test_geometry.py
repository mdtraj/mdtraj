# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

import  os
from mdtraj.testing import get_fn, eq, DocStringFormatTester
import mdtraj.geometry
import numpy as np
from mdtraj.trajectory import load_hdf, load
import mdtraj.trajectory

def test_rg():
    t0 = load(get_fn('frame0.lh5'))
    Rg = mdtraj.geometry.rg.compute_rg(t0)
    Rg0 = np.loadtxt(get_fn("Rg_frame0_ref.dat"))
    eq(Rg, Rg0)

"""
# Compute reference data using MSMBuilder2.
import msmbuilder.geometry.rg
from mdtraj import trajectory

r = trajectory.load("frame0.lh5")
xyz = r.xyz
Rg = msmbuilder.geometry.rg.calculate_rg(xyz)
np.savetxt("Rg_frame0_ref.dat", Rg)
"""

def test_atom_distances():
    t0 = load(get_fn('frame0.lh5'))
    atom_pairs = np.loadtxt(get_fn("atom_pairs.dat"),'int')
    distances = mdtraj.geometry.contact.compute_atom_distances(t0, atom_pairs)
    distances0 = np.loadtxt(get_fn("atom_distances_frame0_ref.dat"))
    eq(distances, distances0)

"""
# Compute reference data using MSMBuilder2.
import msmbuilder.geometry.contact
from mdtraj import trajectory

r = trajectory.load("frame0.lh5")
x = r.xyz
atom_pairs = np.array([[0,1],[0,2],[0,3],[1,2],[2,3],[3,4],[4,5]])
np.savetxt("./atom_pairs.dat",atom_pairs, "%d")
distances = msmbuilder.geometry.contact.atom_distances(x, atom_pairs)
np.savetxt("atom_distances_frame0_ref.dat", distances)
"""

def test_dihedral_indices():    
    traj = load(get_fn('1bpi.pdb'))
    # Manually compare generated indices to known answers.
    phi0_ind = np.array([3, 12, 13, 14]) - 1  # Taken from PDB, so subtract 1
    psi0_ind = np.array([1, 2,  3, 12]) - 1  # Taken from PDB, so subtract 1

    rid, ind = mdtraj.geometry.dihedral._get_indices_phi(traj)
    eq(ind[0], phi0_ind)
    eq(int(rid[0]), 1)

    rid, ind = mdtraj.geometry.dihedral._get_indices_psi(traj)
    eq(ind[0], psi0_ind)
    eq(int(rid[0]), 0)

def test_dihedral_index_offset_generation():    
    traj = load(get_fn('1bpi.pdb'))

    result = np.array([2, 11, 12, 13])  # The atom indices of the first phi angle

    rid1, ind1 = mdtraj.geometry.dihedral._get_indices_phi(traj)
    rid2, ind2 = mdtraj.geometry.dihedral.atom_sequence_finder(traj, ["-C","N","CA","C"])
    rid3, ind3 = mdtraj.geometry.dihedral.atom_sequence_finder(traj, ["-C","N","CA","C"], [-1, 0, 0, 0])
    rid4, ind4 = mdtraj.geometry.dihedral.atom_sequence_finder(traj, ["C" ,"N","CA","C"], [-1, 0, 0, 0])
    eq(rid1, rid2)
    eq(rid1, rid3)
    eq(rid1, rid4)

    eq(ind1, ind2)
    eq(ind1, ind3)
    eq(ind1, ind4)
    eq(ind1[0], result)

def test_dihedral():
    """We compared phi and psi angles from pymol to MDTraj output."""
    traj = load(get_fn('1bpi.pdb'))
    rid, phi = mdtraj.geometry.dihedral.compute_phi(traj)
    phi0 = np.array([-34.50956, -50.869690]) * np.pi / 180.  # Pymol
    eq(phi[0,0:2], phi0, decimal=4)
    eq(int(rid[0]), 1)

    rid, psi = mdtraj.geometry.dihedral.compute_psi(traj)
    psi0 = np.array([134.52554, 144.880173]) * np.pi / 180.  # Pymol
    eq(psi[0,0:2], psi0, decimal=4)
    eq(int(rid[0]), 0)

    rid, chi = mdtraj.geometry.dihedral.compute_chi(traj)
    chi0 = np.array([-43.37841, -18.14592]) * np.pi / 180.  # Pymol
    eq(chi[0,0:2], chi0, decimal=4)
    eq(int(rid[0]), 0)
