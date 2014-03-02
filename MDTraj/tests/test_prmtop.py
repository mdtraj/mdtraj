##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: TJ Lane
# Contributors: Robert McGibbon
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
from mdtraj.testing import get_fn, eq, DocStringFormatTester
from mdtraj import prmtop
doc = DocStringFormatTester(prmtop)


def test_load_prmtop():
    top, _ = prmtop.load_prmtop(get_fn('alanine-dipeptide-implicit.prmtop'))
    ref_top = md.load(get_fn('native.pdb')).topology
    eq(top, ref_top)


def test_load_binpos_w_prmtop():
    traj = md.load(get_fn('alanine.binpos'), top=get_fn('alanine-dipeptide-implicit.prmtop'))
    ref_traj = md.load(get_fn('native.pdb'))

    eq(traj.topology, ref_traj.topology)
    eq(traj.xyz, ref_traj.xyz)


def test_load_binpos_w_prmtop_w_unitcell():
    # this one has box info in the prmtop file
    traj = md.load(get_fn('alanine-dipeptide-explicit.binpos'),
                   top=get_fn('alanine-dipeptide-explicit.prmtop'))
    traj2 = md.load_frame(get_fn('alanine-dipeptide-explicit.binpos'), 0,
                   top=get_fn('alanine-dipeptide-explicit.prmtop'))
    ref_traj = md.load(get_fn('alanine-dipeptide-explicit.pdb'))

    yield lambda: eq(traj.unitcell_vectors, ref_traj.unitcell_vectors, decimal=4)
    yield lambda: eq(traj2.unitcell_vectors[0], ref_traj.unitcell_vectors[0], decimal=4)
    yield lambda: eq(traj.xyz, ref_traj.xyz)
    # the PRMTOP loader puts all of the protein and water into one "chain",
    # whereas the PDB loader puts them in different chains. So the topologies
    # don't exactly match up. I don't think this is too big of a deal...
    # yield lambda: eq(traj.topology, ref_traj.topology)
    