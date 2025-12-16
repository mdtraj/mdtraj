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
# License along with  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import math
import tempfile

import numpy as np

from mdtraj import element, load
from mdtraj.core.topology import Topology
from mdtraj.core.trajectory import Trajectory
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

            traj1.save(file.name)

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


def create_test_topology():
    """Create a test topology with 3 chains, including a disulfide bond and a non-standard residue."""
    topology = Topology()
    # Chain 1
    chain1 = topology.add_chain("A")
    res1a = topology.add_residue("GLY", chain1, 1)
    n1 = topology.add_atom("N", element.nitrogen, res1a)
    ca1 = topology.add_atom("CA", element.carbon, res1a)
    c1 = topology.add_atom("C", element.carbon, res1a)
    topology.add_bond(n1, ca1)
    topology.add_bond(ca1, c1)
    res2a = topology.add_residue("ALA", chain1, 2)
    n2 = topology.add_atom("N", element.nitrogen, res2a)
    ca2 = topology.add_atom("CA", element.carbon, res2a)
    c2 = topology.add_atom("C", element.carbon, res2a)
    topology.add_bond(c1, n2)
    topology.add_bond(n2, ca2)
    topology.add_bond(ca2, c2)
    res3a = topology.add_residue("CYS", chain1, 3)
    n3 = topology.add_atom("N", element.nitrogen, res3a)
    ca3 = topology.add_atom("CA", element.carbon, res3a)
    cb3 = topology.add_atom("CB", element.carbon, res3a)
    c3 = topology.add_atom("C", element.carbon, res3a)
    sg1 = topology.add_atom("SG", element.sulfur, res3a)
    topology.add_bond(c2, n3)
    topology.add_bond(n3, ca3)
    topology.add_bond(ca3, c3)
    topology.add_bond(ca3, cb3)
    topology.add_bond(cb3, sg1)
    # Chain 2
    chain2 = topology.add_chain("B")
    res1b = topology.add_residue("GLY", chain2, 1)
    n4 = topology.add_atom("N", element.nitrogen, res1b)
    n4h1 = topology.add_atom("H", element.nitrogen, res1b)
    ca4 = topology.add_atom("CA", element.carbon, res1b)
    cah1 = topology.add_atom("HA2", element.hydrogen, res1b)
    cah2 = topology.add_atom("HA3", element.hydrogen, res1b)
    c4 = topology.add_atom("C", element.carbon, res1b)
    o4 = topology.add_atom("O", element.oxygen, res1b)
    topology.add_bond(n4, ca4)
    topology.add_bond(ca4, c4)
    topology.add_bond(ca4, cah1)
    topology.add_bond(ca4, cah2)
    topology.add_bond(c4, o4)
    topology.add_bond(n4, n4h1)
    res2b = topology.add_residue("ALA", chain2, 2)
    n5 = topology.add_atom("N", element.nitrogen, res2b)
    ca5 = topology.add_atom("CA", element.carbon, res2b)
    cb5 = topology.add_atom("CB", element.carbon, res2b)
    c5 = topology.add_atom("C", element.carbon, res2b)
    o5 = topology.add_atom("O", element.oxygen, res2b)
    h1 = topology.add_atom("H", element.hydrogen, res2b)
    h2 = topology.add_atom("HA", element.hydrogen, res2b)
    h3 = topology.add_atom("HB1", element.hydrogen, res2b)
    h4 = topology.add_atom("HB2", element.hydrogen, res2b)
    h5 = topology.add_atom("HB3", element.hydrogen, res2b)
    topology.add_bond(c4, n5)
    topology.add_bond(n5, h1)
    topology.add_bond(ca5, h2)
    topology.add_bond(ca5, cb5)
    topology.add_bond(cb5, h3)
    topology.add_bond(cb5, h4)
    topology.add_bond(cb5, h5)
    topology.add_bond(ca5, c5)
    topology.add_bond(ca5, n5)
    topology.add_bond(c5, o5)
    res3b = topology.add_residue("CYS", chain2, 3)
    n6 = topology.add_atom("N", element.nitrogen, res3b)
    ca6 = topology.add_atom("CA", element.carbon, res3b)
    cb6 = topology.add_atom("CB", element.carbon, res3b)
    c6 = topology.add_atom("C", element.carbon, res3b)
    sg2 = topology.add_atom("SG", element.sulfur, res3b)
    o61 = topology.add_atom("O", element.oxygen, res3b)
    o62 = topology.add_atom("OXT", element.oxygen, res3b)
    topology.add_bond(c5, n6)
    topology.add_bond(n6, ca6)
    topology.add_bond(ca6, c6)
    topology.add_bond(ca6, cb6)
    topology.add_bond(cb6, sg2)
    topology.add_bond(c6, o61)
    topology.add_bond(c6, o62)
    topology.add_bond(sg1, sg2)
    # Chain 3
    chain3 = topology.add_chain("C")
    res1c = topology.add_residue("2HA", chain3, 1)
    cc1 = topology.add_atom("C1", element.carbon, res1c)
    cc2 = topology.add_atom("C2", element.carbon, res1c)
    cc3 = topology.add_atom("C3", element.carbon, res1c)
    ch11 = topology.add_atom("1H1C", element.hydrogen, res1c)
    ch12 = topology.add_atom("2H1C", element.hydrogen, res1c)
    ch31 = topology.add_atom("1H3C", element.hydrogen, res1c)
    ch32 = topology.add_atom("2H3C", element.hydrogen, res1c)
    ho1 = topology.add_atom("H1", element.hydrogen, res1c)
    ho3 = topology.add_atom("H3", element.hydrogen, res1c)
    co1 = topology.add_atom("O1", element.oxygen, res1c)
    co2 = topology.add_atom("O2", element.oxygen, res1c)
    co3 = topology.add_atom("O3", element.oxygen, res1c)
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
    return topology


def assert_traj_fields_equal(traj1, traj2):
    assert np.allclose(traj1.xyz, traj2.xyz, atol=1e-5), "Positions do not match after write/read"
    assert traj1.n_atoms == traj2.n_atoms, "Number of atoms does not match"
    assert traj1.n_residues == traj2.n_residues, "Number of residues does not match"
    assert traj1.n_chains == traj2.n_chains, "Number of chains does not match"
    assert np.allclose(traj1.unitcell_lengths, traj2.unitcell_lengths, atol=1e-5), "Unit cell lengths do not match"
    assert np.allclose(traj1.unitcell_angles, traj2.unitcell_angles, atol=1e-5), "Unit cell angles do not match"


def bond_set(top):
    return set(tuple(sorted((a1.index, a2.index))) for a1, a2 in top.bonds)


def assert_bonds_equal(top1, top2):
    bonds1 = bond_set(top1)
    bonds2 = bond_set(top2)
    assert bonds1 == bonds2, f"Bonds do not match after write/read: {bonds1} != {bonds2}"


def test_write_and_read_pdbx(tmp_path):
    """Test writing a PDBx file and reading it back,
    checking all fields including bonds, with two chains and a disulfide bond."""
    topology = create_test_topology()
    n_atoms = topology.n_atoms
    xyz = np.zeros((1, n_atoms, 3), dtype=np.float32)
    for i in range(n_atoms):
        xyz[0, i, :] = [i * 0.1, 0.0, 0.0]
    unitcell_lengths = np.array([[2.0, 2.0, 2.0]], dtype=np.float32)
    unitcell_angles = np.array([[90.0, 90.0, 90.0]], dtype=np.float32)
    traj = Trajectory(xyz=xyz, topology=topology, unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles)
    out_path = tmp_path / "test_write.pdbx"
    traj.save(str(out_path))
    loaded = load(str(out_path))
    assert_traj_fields_equal(traj, loaded)
    assert_bonds_equal(traj.topology, loaded.topology)


def test_write_and_read_pdbx_4py5(tmp_path, get_fn):
    """Test writing and reading 4py5.cif, checking all fields including bonds."""
    cif_path = get_fn("4py5.cif")
    traj = load(cif_path)
    out_path = tmp_path / "test_write_4py5.pdbx"
    traj.save(str(out_path))
    loaded = load(str(out_path))
    assert_traj_fields_equal(traj, loaded)
    assert len(bond_set(traj.topology)) == 2928
    assert_bonds_equal(traj.topology, loaded.topology)


def test_write_and_read_pdbx_from_dcd_with_bonds_and_keepids(tmp_path, get_fn):
    """Test writing CIF from DCD+PDB topology with bonds, ensuring bonds are preserved.

    This test addresses the issue where bonds were lost when converting DCD trajectories
    with PDB topologies to CIF format. The problem occurred because:
    1. Chain and residue objects in bonds were different instances than those in topology.chains/residues
    2. Using objects as dictionary keys failed due to Python's object identity hashing

    The fix uses .index properties (unique integers) as dictionary keys instead of objects.
    """
    # Load DCD with PDB topology that has bonds
    dcd_path = get_fn("ala_dipeptide.dcd")
    pdb_path = get_fn("ala_dipeptide.pdb")
    traj = load(dcd_path, top=pdb_path)

    # Verify we have bonds in the original topology
    orig_bonds = len(list(traj.topology.bonds))
    assert orig_bonds > 0, "Original topology has no bonds to test preservation"

    # Save first frame to CIF with keepIds=False (default - auto-generate chain IDs)
    cif_path_no_keepids = tmp_path / "test_dcd_to_cif_no_keepids.cif"
    traj[0].save_cif(str(cif_path_no_keepids), keepIds=False)

    # Load back and verify bonds are preserved
    loaded_no_keepids = load(str(cif_path_no_keepids))
    loaded_bonds_no_keepids = len(list(loaded_no_keepids.topology.bonds))
    assert loaded_bonds_no_keepids == orig_bonds, (
        f"Bonds lost with keepIds=False: {orig_bonds} -> {loaded_bonds_no_keepids}"
    )
    assert_bonds_equal(traj[0].topology, loaded_no_keepids.topology)

    # Save with keepIds=True (preserve original chain IDs from PDB)
    cif_path_keepids = tmp_path / "test_dcd_to_cif_keepids.cif"
    traj[0].save_cif(str(cif_path_keepids), keepIds=True)

    # Load back and verify bonds are preserved
    loaded_keepids = load(str(cif_path_keepids))
    loaded_bonds_keepids = len(list(loaded_keepids.topology.bonds))
    assert loaded_bonds_keepids == orig_bonds, f"Bonds lost with keepIds=True: {orig_bonds} -> {loaded_bonds_keepids}"
    assert_bonds_equal(traj[0].topology, loaded_keepids.topology)


def test_residue_selection_from_cif(get_fn):
    """Tests that residues are correctly parsed as integers in CIF files"""

    # Load the trajectory from the test CIF file using the fixture
    cif_file = get_fn("8ddg.cif")
    traj = load(cif_file)
    top = traj.topology

    # Try selecting by residue number
    residue_atoms = top.select("residue 67")  # ResSeq based index
    resid_atoms = top.select("resid 0")  # 0-based index

    assert residue_atoms.size > 0, "Expected atoms in residue 67 from CIF"
    assert len(residue_atoms) == len(resid_atoms), "Expected same number of atoms in residue 67 and resid 0"
    assert np.array_equal(residue_atoms, resid_atoms), "Expected atoms in residue 67 and resid 0 to be the same"
    assert np.array_equal(
        top.select("residue 68"),
        top.select("resid 1"),
    ), "Expected atoms in residue 68 and resid 1 to be the same"
