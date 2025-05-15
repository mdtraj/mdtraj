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

import math
import tempfile

import numpy as np

from mdtraj import load
from mdtraj.formats import PDBxTrajectoryFile
from mdtraj.testing import eq


def test_convert(get_fn):
    for filename in ["2EQQ.pdb", "4OH9.pdb"]:
        # Load a PDB file.

        traj1 = load(get_fn(filename))
        with tempfile.NamedTemporaryFile(
            suffix=".pdbx",
            mode="w",
            delete=False,
        ) as file:
            # Save it in PDBx/mmCIF format.

            pdbx1 = PDBxTrajectoryFile(file.name, mode="w")
            for i in range(traj1.n_frames):
                pdbx1.write(
                    traj1.xyz[i],
                    traj1.topology,
                    traj1.unitcell_lengths,
                    traj1.unitcell_angles,
                )
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
            assert eq(traj3.n_frames, math.ceil(traj1.n_frames / 2))
            assert eq(traj3.n_atoms, 10)
            assert eq(traj1.xyz[::2, 10:20], traj3.xyz)
            atoms1 = list(traj1.topology.atoms)[10:20]
            for a1, a2 in zip(atoms1, traj3.topology.atoms):
                assert eq(a1.name, a2.name)


def test_load_pdbx_file(get_fn):
    """Test loading a PDBx/mmCIF file and verifying topology details."""
    for filename in ["8ddg.cif", "8ddg.cif.gz"]:
        traj = load(get_fn(filename))

        # Verify the number of atoms, residues, chains, and bonds
        assert len([a for a in traj.topology.atoms]) == 64, "Unexpected number of atoms"
        assert len([r for r in traj.topology.residues]) == 3, "Unexpected number of residues"
        assert len([c for c in traj.topology.chains]) == 1, "Unexpected number of chains"
        assert len([b for b in traj.topology.bonds]) == 66, "Unexpected number of bonds"

        # Verify atom names and residue names
        atom_names = [a.name for a in traj.topology.atoms]
        assert "CA" in atom_names, "Expected atom 'CA' not found"
        assert "OXT" in atom_names, "Expected atom 'OXT' not found"
        assert "HE2" in atom_names, "Expected atom 'HE2' not found"
        assert "HXT" not in atom_names, "Unexpected atom 'HXT' found"

        residue_names = [r.name for r in traj.topology.residues]
        assert "TYR" in residue_names, "Expected residue 'TYR' not found"
        assert "PHE" in residue_names, "Expected residue 'PHE' not found"
        assert "ALA" not in residue_names, "Unexpected residue 'ALA' found"

        # Verify chain IDs
        chain_ids = [c.chain_id for c in traj.topology.chains]
        assert "A" in chain_ids, "Expected chain ID 'A' not found"

        # Verify unit cell lengths and angles
        assert traj.unitcell_lengths is not None, "Unit cell lengths not found"
        assert traj.unitcell_angles is not None, "Unit cell angles not found"

        # Expected values for unit cell lengths (in nanometers) and angles (in degrees)
        expected_lengths = [
            2.3140,
            0.4840,
            1.9790,
        ]  # Converted from Ångstroms to nanometers
        expected_angles = [90.000, 107.048, 90.000]

        # Use numpy's allclose for comparison with a tolerance
        assert np.allclose(
            traj.unitcell_lengths[0],
            expected_lengths,
            atol=1e-3,
        ), "Unit cell lengths do not match"
        assert np.allclose(
            traj.unitcell_angles[0],
            expected_angles,
            atol=1e-3,
        ), "Unit cell angles do not match"
