##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2022 Stanford University and the Authors
#
# Authors: Peter Eastman
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

from mdtraj.formats import PDBxTrajectoryFile
from mdtraj.testing import eq
from mdtraj import load
import math
import tempfile
import pytest

try:
    import openmm
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

# special pytest global to mark all tests in this module
pytestmark = pytest.mark.skipif(not HAVE_OPENMM, reason='test_pdbx.py needs OpenMM.')

def test_convert(get_fn):
    for filename in ['2EQQ.pdb', '4OH9.pdb']:
        # Load a PDB file.

        traj1 = load(get_fn(filename))
        with tempfile.NamedTemporaryFile(suffix='.pdbx', mode='w', delete=False) as file:
            # Save it in PDBx/mmCIF format.

            pdbx1 = PDBxTrajectoryFile(file.name, mode='w')
            for i in range(traj1.n_frames):
                pdbx1.write(traj1.xyz[i], traj1.topology, traj1.unitcell_lengths, traj1.unitcell_angles)
            pdbx1.close()

            # Load the PDBx/mmCIF file and make the result is identical.

            traj2 = load(file.name)
            assert eq(traj1.n_frames, traj2.n_frames)
            assert eq(traj1.n_atoms, traj2.n_atoms)
            assert eq(traj1.xyz, traj2.xyz)
            assert eq(traj1.unitcell_lengths, traj2.unitcell_lengths)
            assert eq(traj1.unitcell_angles, traj2.unitcell_angles)
            for a1, a2 in zip(traj1.topology.atoms, traj2.topology.atoms):
                assert eq(a1, a2)

            # Try loading just a subset of the frames and atoms.

            traj3 = load(file.name, atom_indices=range(10, 20), stride=2)
            assert eq(traj3.n_frames, math.ceil(traj1.n_frames/2))
            assert eq(traj3.n_atoms, 10)
            assert eq(traj1.xyz[::2, 10:20], traj3.xyz)
            atoms1 = list(traj1.topology.atoms)[10:20]
            for a1, a2 in zip(atoms1, traj3.topology.atoms):
                assert eq(a1.name, a2.name)
