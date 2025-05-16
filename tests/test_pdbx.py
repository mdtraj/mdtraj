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

import pytest

from mdtraj import load
from mdtraj.formats import PDBxTrajectoryFile
from mdtraj.testing import eq
import numpy as np


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
        assert (
            len([r for r in traj.topology.residues]) == 3
        ), "Unexpected number of residues"
        assert (
            len([c for c in traj.topology.chains]) == 1
        ), "Unexpected number of chains"
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
            traj.unitcell_lengths[0], expected_lengths, atol=1e-3
        ), "Unit cell lengths do not match"
        assert np.allclose(
            traj.unitcell_angles[0], expected_angles, atol=1e-3
        ), "Unit cell angles do not match"


def test_write_and_read_pdbx(tmp_path):
    """Test writing a PDBx file and reading it back, checking all fields including bonds, with two chains and a disulfide bond."""
    import mdtraj as md
    import numpy as np
    from mdtraj import Trajectory
    from mdtraj.formats import PDBxTrajectoryFile

    # Create topology: 2 chains, 3 residues each
    topology = md.Topology()
    # Chain 1: backbone only (N, CA, C)
    chain1 = topology.add_chain("A")
    res1a = topology.add_residue("GLY", chain1, 1)
    n1 = topology.add_atom("N", md.element.nitrogen, res1a) # 0
    ca1 = topology.add_atom("CA", md.element.carbon, res1a) # 1
    c1 = topology.add_atom("C", md.element.carbon, res1a) # 2
    topology.add_bond(n1, ca1)
    topology.add_bond(ca1, c1)
    res2a = topology.add_residue("ALA", chain1, 2)
    n2 = topology.add_atom("N", md.element.nitrogen, res2a) # 3
    ca2 = topology.add_atom("CA", md.element.carbon, res2a) # 4    
    c2 = topology.add_atom("C", md.element.carbon, res2a) # 5
    topology.add_bond(c1, n2) # peptide bond
    topology.add_bond(n2, ca2)
    topology.add_bond(ca2, c2)
    res3a = topology.add_residue("CYS", chain1, 3)
    n3 = topology.add_atom("N", md.element.nitrogen, res3a) # 6
    ca3 = topology.add_atom("CA", md.element.carbon, res3a) # 7
    cb3 = topology.add_atom("CB", md.element.carbon, res3a) # 8
    c3 = topology.add_atom("C", md.element.carbon, res3a) # 9
    sg1 = topology.add_atom("SG", md.element.sulfur, res3a) # 10
    topology.add_bond(c2, n3) # peptide bond
    topology.add_bond(n3, ca3)
    topology.add_bond(ca3, c3)
    topology.add_bond(ca3, cb3)
    topology.add_bond(cb3, sg1)

    # Chain 2: all heavy atoms in first residue, all hydrogens in second, CYS with SG in third
    chain2 = topology.add_chain("B")
    # Residue 1: all heavy atoms (GLY)
    res1b = topology.add_residue("GLY", chain2, 1)
    n4 = topology.add_atom("N", md.element.nitrogen, res1b) # 11
    n4h1 = topology.add_atom("H", md.element.nitrogen, res1b) # 12
    ca4 = topology.add_atom("CA", md.element.carbon, res1b) # 13
    cah1 = topology.add_atom("HA2", md.element.hydrogen, res1b) # 14
    cah2 = topology.add_atom("HA3", md.element.hydrogen, res1b) # 14
    c4 = topology.add_atom("C", md.element.carbon, res1b) # 15
    o4 = topology.add_atom("O", md.element.oxygen, res1b) # 16
    topology.add_bond(n4, ca4)
    topology.add_bond(ca4, c4)
    topology.add_bond(ca4, cah1)
    topology.add_bond(ca4, cah2)
    topology.add_bond(c4, o4)
    topology.add_bond(n4, n4h1)
    
    
    # Residue 2: all hydrogens (ALA)
    res2b = topology.add_residue("ALA", chain2, 2)
    n5 = topology.add_atom("N", md.element.nitrogen, res2b) # 17
    ca5 = topology.add_atom("CA", md.element.carbon, res2b) # 18
    cb5 = topology.add_atom("CB", md.element.carbon, res2b) # 19
    c5 = topology.add_atom("C", md.element.carbon, res2b) # 20
    o5 = topology.add_atom("O", md.element.oxygen, res2b) # 21
    h1 = topology.add_atom("H", md.element.hydrogen, res2b) # 22
    h2 = topology.add_atom("HA", md.element.hydrogen, res2b) # 23
    h3 = topology.add_atom("HB1", md.element.hydrogen, res2b) # 24
    h4 = topology.add_atom("HB2", md.element.hydrogen, res2b) # 25
    h5 = topology.add_atom("HB3", md.element.hydrogen, res2b) # 26
    topology.add_bond(c4, n5) # peptide bond
    topology.add_bond(n5, h1)
    topology.add_bond(ca5, h2)
    topology.add_bond(ca5, cb5)
    topology.add_bond(cb5, h3)
    topology.add_bond(cb5, h4)
    topology.add_bond(cb5, h5)
    topology.add_bond(ca5, c5)
    topology.add_bond(ca5, n5)
    topology.add_bond(c5, o5)


    # Residue 3: CYS with SG
    res3b = topology.add_residue("CYS", chain2, 3)
    n6 = topology.add_atom("N", md.element.nitrogen, res3b) # 27
    ca6 = topology.add_atom("CA", md.element.carbon, res3b) # 28
    cb6 = topology.add_atom("CB", md.element.carbon, res3b) # 29
    c6 = topology.add_atom("C", md.element.carbon, res3b) # 30
    sg2 = topology.add_atom("SG", md.element.sulfur, res3b) # 31
    o61 = topology.add_atom("O", md.element.oxygen, res3b) # 32
    o62 = topology.add_atom("OXT", md.element.oxygen, res3b) # 33
    oh6 = topology.add_atom("OH", md.element.hydrogen, res3b) # 34
    topology.add_bond(c5, n6) # peptide bond
    topology.add_bond(n6, ca6)
    topology.add_bond(ca6, c6)
    topology.add_bond(ca6, cb6)
    topology.add_bond(cb6, sg2)
    topology.add_bond(c6, o61)
    topology.add_bond(c6, o62)



    # Disulfide bond between the two CYS SG atoms
    topology.add_bond(sg1, sg2) #Disulfide bond

    # Chain C: A new residue
    chain3 = topology.add_chain("C")
    res1c = topology.add_residue("2HA", chain3, 1)
    cc1 = topology.add_atom("C1", md.element.carbon, res1c) # 32
    cc2 = topology.add_atom("C2", md.element.carbon, res1c) # 33
    cc3 = topology.add_atom("C3", md.element.carbon, res1c) # 34
    ch11 = topology.add_atom("1H1C", md.element.hydrogen, res1c) # 35
    ch12 = topology.add_atom("2H1C", md.element.hydrogen, res1c) # 36
    ch31 = topology.add_atom("1H3C", md.element.hydrogen, res1c) # 37
    ch32 = topology.add_atom("2H3C", md.element.hydrogen, res1c) # 38
    ho1 = topology.add_atom("H1", md.element.hydrogen, res1c) # 39
    ho3 = topology.add_atom("H3", md.element.hydrogen, res1c) # 40
    co1 = topology.add_atom("O1", md.element.oxygen, res1c) # 41
    co2 = topology.add_atom("O2", md.element.oxygen, res1c) # 42
    co3 = topology.add_atom("O3", md.element.oxygen, res1c) # 43
    topology.add_bond(cc1, cc2)
    topology.add_bond(cc2, cc3)
    topology.add_bond(cc1, ch11)
    topology.add_bond(cc1, ch12)
    topology.add_bond(cc2, ch31)
    topology.add_bond(cc2, ch32)
    topology.add_bond(cc1, ho1)
    topology.add_bond(cc3, ho3)
    topology.add_bond(cc1, co1)
    topology.add_bond(cc2, co2)
    topology.add_bond(cc3, co3)

    # Create coordinates for all atoms (1 frame)
    n_atoms = topology.n_atoms
    xyz = np.zeros((1, n_atoms, 3), dtype=np.float32)
    for i in range(n_atoms):
        xyz[0, i, :] = [i * 0.1, 0.0, 0.0]
    unitcell_lengths = np.array([[2.0, 2.0, 2.0]], dtype=np.float32)
    unitcell_angles = np.array([[90.0, 90.0, 90.0]], dtype=np.float32)
    traj = Trajectory(xyz=xyz, topology=topology, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)

    # Write to a temporary PDBx file
    from pathlib import Path
    tmp_path= Path('.')
    out_path = tmp_path / "test_write.pdbx"
    with PDBxTrajectoryFile(str(out_path), mode="w") as f:
        f.write(traj.xyz[0], traj.topology, traj.unitcell_lengths, traj.unitcell_angles)

    # Uncomment the following line to print the output file path for inspection
    # print(f"PDBx file written to: {out_path}")

    # Read it back
    loaded = md.load(str(out_path))

    # Check positions
    assert np.allclose(loaded.xyz, traj.xyz, atol=1e-5), "Positions do not match after write/read"
    # Check topology: atoms, residues, chains
    assert loaded.n_atoms == traj.n_atoms, "Number of atoms does not match"
    assert loaded.n_residues == traj.n_residues, "Number of residues does not match"
    assert loaded.n_chains == traj.n_chains, "Number of chains does not match"
    # Check unit cell
    assert np.allclose(loaded.unitcell_lengths, traj.unitcell_lengths, atol=1e-5), "Unit cell lengths do not match"
    assert np.allclose(loaded.unitcell_angles, traj.unitcell_angles, atol=1e-5), "Unit cell angles do not match"

    # Helper to classify bond type
    def classify_bond(a1, a2):
        res1, res2 = a1.residue, a2.residue
        # Disulfide bond: CYS SG-SG between different residues/chains
        if (
            res1.name == "CYS" and res2.name == "CYS" and a1.name == "SG" and a2.name == "SG" and res1 != res2
        ):
            return "disulfide"
        # Peptide bond: C-N between adjacent residues in the same chain
        if (
            a1.name == "C" and a2.name == "N" and res1.chain == res2.chain and abs(int(res1.resSeq) - int(res2.resSeq)) == 1
        ) or (
            a2.name == "C" and a1.name == "N" and res1.chain == res2.chain and abs(int(res1.resSeq) - int(res2.resSeq)) == 1
        ):
            return "peptide"
        # Common residue bonds: backbone N-CA, CA-C, C-O, CA-CB, etc. (same residue, common names)
        common_pairs = {frozenset(["N", "CA"]), frozenset(["CA", "C"]), frozenset(["C", "O"]), frozenset(["CA", "CB"]), frozenset(["CA", "SG"])}
        if res1 == res2 and frozenset([a1.name, a2.name]) in common_pairs:
            return "common"
        # Uncommon residue bonds: same residue, but not in common_pairs
        if res1 == res2:
            return "uncommon"
        # Otherwise, unknown
        return "other"

    def bonds_by_type(top):
        types = {"common": set(), "uncommon": set(), "peptide": set(), "disulfide": set(), "other": set()}
        for a1, a2 in top.bonds:
            btype = classify_bond(a1, a2)
            idx_pair = tuple(sorted((a1.index, a2.index)))
            types[btype].add(idx_pair)
        return types

    bonds_written = bonds_by_type(traj.topology)
    bonds_read = bonds_by_type(loaded.topology)

    debug_msgs = []
    all_passed = True
    for btype in ["common", "uncommon", "peptide", "disulfide"]:
        written = bonds_written[btype]
        read = bonds_read[btype]
        missing_in_read = written - read
        extra_in_read = read - written
        if missing_in_read or extra_in_read:
            all_passed = False
            if missing_in_read:
                for idx_pair in missing_in_read:
                    a1 = traj.topology.atom(idx_pair[0])
                    a2 = traj.topology.atom(idx_pair[1])
                    debug_msgs.append(
                        f"[{btype}] Written but not read ({a1.index}-{a2.index}): Bond between {a1.name} ({a1.residue.chain.chain_id} {a1.residue.name} {a1.residue.resSeq}) and {a2.name} ({a1.residue.chain.chain_id} {a2.residue.name} {a2.residue.resSeq}) [type: {btype}]"
                    )
            if extra_in_read:
                for idx_pair in extra_in_read:
                    a1 = loaded.topology.atom(idx_pair[0])
                    a2 = loaded.topology.atom(idx_pair[1])
                    debug_msgs.append(
                        f"[{btype}] Read but not written ({a1.index}-{a2.index}): Bond between {a1.name} ({a1.residue.chain.chain_id} {a1.residue.name} {a1.residue.resSeq}) and {a2.name} ({a1.residue.chain.chain_id} {a2.residue.name} {a2.residue.resSeq}) [type: {btype}]"
                    )
    if not all_passed:
        debug_str = "\n".join(debug_msgs)
        raise AssertionError(f"Bond mismatch after write/read. Details:\n{debug_str}")
