##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2025 Stanford University and the Authors
#
# Authors: Lee-Ping Wang
# Contributors: Robert McGibbon, Jeremy M. G. Leung
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


import numpy as np

import mdtraj as md
from mdtraj.formats import ArcTrajectoryFile, PDBTrajectoryFile
from mdtraj.formats.arc import load_arc
from mdtraj.testing import eq


def test_read_0(get_fn):
    """Check if .arc positions are read correctly"""
    with ArcTrajectoryFile(get_fn("4waters.arc")) as f:
        xyz, leng, ang = f.read()
    with PDBTrajectoryFile(get_fn("4waters.pdb")) as f:
        xyz2 = f.positions
    eq(xyz, xyz2, decimal=3)


def test_read_arctraj(get_fn):
    traj = md.load(get_fn("nitrogen.arc"), top=get_fn("nitrogen.pdb"))
    owntop = md.load(get_fn("nitrogen.arc"))
    eq(traj.xyz, owntop.xyz)


def test_read_arctraj_top(get_fn):
    """Check if anything other than residues differ between loading
    a ``.arc`` file with two different topologies"""
    pdbtop = load_arc(get_fn("4waters.arc"), top=get_fn("4waters.pdb"))
    arctop = load_arc(get_fn("4waters.arc"))

    eq(pdbtop.n_atoms, arctop.n_atoms)
    eq(pdbtop.xyz, arctop.xyz)

    # Only the residues should differ
    assert pdbtop.n_residues == 4
    assert arctop.n_residues == 1

    for ipdb, iarc in zip(pdbtop.top._atoms, arctop.top._atoms):
        # check if atoms are in the same order with same element identity
        eq(ipdb.name, iarc.name)
        eq(ipdb.element, iarc.element)

    for ipdb, iarc in zip(pdbtop.top._bonds, arctop.top._bonds):
        # Check if connectivity is the same, based on index, which is
        # implicitly checked in the for loop above
        eq(ipdb.atom1.index, iarc.atom1.index)
        eq(ipdb.atom2.index, iarc.atom2.index)


def test_read_pbc(get_fn):
    with ArcTrajectoryFile(get_fn("thf200.arc")) as f:
        xyz, leng, ang = f.read()
    eq(xyz[0, 0], np.array([21.231166, 32.549819, 8.278454]), decimal=3)
    eq(xyz[0, -1], np.array([29.013880, 6.255754, 44.074519]), decimal=3)
    eq(leng, np.array([[45.119376, 45.119376, 45.119376]]), decimal=3)
    eq(ang, np.array([[90.0, 90.0, 90.0]]), decimal=3)
