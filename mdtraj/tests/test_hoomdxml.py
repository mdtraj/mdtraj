##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Tim Moore
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

from mdtraj.testing import get_fn, eq
from mdtraj import load


def test_load():
    t = load(get_fn('test_hoomdxml.dcd'), top=get_fn('well-mixed.hoomdxml'))
    top = load(get_fn('well-mixed.hoomdxml')).topology
    eq(t.topology, top)
    eq(top.n_atoms, 6600)
    eq(top.n_bonds, 3200)
    eq(top.n_chains, 4000)
    eq(top.n_chains, top.n_residues)  # hoomdxml assumes each chain is 1 residue
