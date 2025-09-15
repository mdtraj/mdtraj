##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: TJ Lane
# Contributors: Robert McGibbon, Jason Swails
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
from mdtraj.formats import prmtop
from mdtraj.testing import eq


def test_load_prmtop(get_fn):
    top = prmtop.load_prmtop(get_fn("alanine-dipeptide-implicit.prmtop"))
    ref_top = md.load(get_fn("native.pdb")).topology
    eq(top, ref_top)


def test_load_chamber_prmtop(get_fn):
    top = prmtop.load_prmtop(get_fn("ala3_chamber.parm7"))
    eq(top.n_atoms, 33)
    eq(top.n_residues, 3)
    eq(top.n_bonds, 32)
    eq(top.n_chains, 1)
