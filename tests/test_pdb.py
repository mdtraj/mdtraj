##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Jason Swails, Matthew Harrigan
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


import os
import re
import tempfile

import numpy as np
import pytest
from conftest import flaky_pdb_dl

from mdtraj import Topology, load, load_pdb
from mdtraj.core.topology import float_to_bond_type
from mdtraj.formats.pdb import pdbstructure
from mdtraj.formats.pdb.pdbstructure import PdbStructure
from mdtraj.testing import eq

fd, temp = tempfile.mkstemp(suffix=".pdb")
os.close(fd)


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.unlink(temp)


def test_pdbread(get_fn):
    load(get_fn("native.pdb"))


def test_pdbread_with_input_top(get_fn):
    pdb = get_fn("native.pdb")
    p_1 = load(pdb)

    p_2 = load(pdb, top=pdb)

    eq(p_1.xyz, p_2.xyz)
    eq(p_1.topology, p_2.topology)


def test_pdbwrite(get_fn):
    pdb = get_fn("native.pdb")
    p = load(pdb)
    p.save(temp)

    r = load(temp)
    eq(p.xyz, r.xyz)


def test_load_multiframe(get_fn):
    with open(get_fn("multiframe.pdb")) as f:
        pdb = PdbStructure(f)
        assert eq(len(pdb.models), 2)
        assert eq(len(pdb.models[0].chains), 1)
        assert eq(len(pdb.models[0].chains[0].residues), 3)
        assert eq(len(list(pdb.models[0].iter_atoms())), 22)

        assert eq(len(pdb.models[1].chains), 1)
        assert eq(len(pdb.models[1].chains[0].residues), 3)
        assert eq(len(list(pdb.models[1].iter_atoms())), 22)

    t = load(get_fn("multiframe.pdb"))
    assert eq(t.n_frames, 2)
    assert eq(t.n_atoms, 22)
    assert eq(t.xyz[0], t.xyz[1])


def test_4ZUO(get_fn):
    t = load(get_fn("4ZUO.pdb"))
    eq(t.n_frames, 1)
    eq(t.n_atoms, 6200)

    # this is a random line from the file
    # ATOM   1609  O   GLY A 201     -25.423  13.774 -25.234  1.00  8.92           O

    atom = list(t.top.atoms)[1525]
    eq(atom.element.symbol, "O")
    eq(atom.residue.name, "GLY")
    eq(atom.residue.chain.chain_id, "A")
    eq(atom.index, 1525)
    eq(t.xyz[0, 1525], np.array([-25.423, 13.774, -25.234]) / 10)  # converting to NM

    # this is atom 1577 in the PDB
    # CONECT 1577 5518
    # ATOM   1577  O   HIS A 197     -18.521   9.724 -32.805  1.00  8.81           O
    # HETATM 5518  K     K A 402     -19.609  10.871 -35.067  1.00  9.11           K
    atom = list(t.top.atoms)[1493]
    eq(atom.name, "O")
    eq(atom.residue.name, "HIS")
    eq(
        [(a1.index, a2.index) for a1, a2 in t.top.bonds if a1.index == 1493 or a2.index == 1493],
        [(1492, 1493), (1493, 5129)],
    )

    # that first bond is from a conect record


def test_2EQQ_0(get_fn):
    # this is an nmr structure with 20 models
    t = load(get_fn("2EQQ.pdb"))
    assert eq(t.n_frames, 20)
    assert eq(t.n_atoms, 423)
    assert eq(len(list(t.top.residues)), 28)


def test_1vii_solvated_with_ligand(get_fn):
    traj = load(get_fn("1vii_sustiva_water.pdb"))
    eq(len(list(traj.top.bonds)), 5156)
    eq(len([bond for bond in traj.top.bonds if bond[0].residue.name == "LIG"]), 32)
    traj.save(temp)
    traj = load(temp)
    eq(len(list(traj.top.bonds)), 5156)
    eq(len([bond for bond in traj.top.bonds if bond[0].residue.name == "LIG"]), 32)


def test_write_large(get_fn):
    traj = load(get_fn("native.pdb"))
    traj.xyz.fill(123456789)
    with pytest.raises(ValueError):
        traj.save(temp)


def test_write_large_2(get_fn):
    traj = load(get_fn("native.pdb"))
    traj.xyz.fill(-123456789)
    with pytest.raises(ValueError):
        traj.save(temp)


def test_pdbstructure_0():
    pdb_lines = [
        "ATOM    188  N   CYS A  42      40.714  -5.292  12.123  1.00 11.29           N  ",
        "ATOM    189  CA  CYS A  42      39.736  -5.883  12.911  1.00 10.01           C  ",
        "ATOM    190  C   CYS A  42      40.339  -6.654  14.087  1.00 22.28           C  ",
        "ATOM    191  O   CYS A  42      41.181  -7.530  13.859  1.00 13.70           O  ",
        "ATOM    192  CB  CYS A  42      38.949  -6.825  12.002  1.00  9.67           C  ",
        "ATOM    193  SG  CYS A  42      37.557  -7.514  12.922  1.00 20.12           S  ",
    ]

    res = pdbstructure.Residue("CYS", 42)
    for line in pdb_lines:
        res._add_atom(pdbstructure.Atom(line))
    for i, atom in enumerate(res):
        eq(pdb_lines[i], str(atom))


def test_pdbstructure_1():
    pdb_lines = [
        "ATOM    188  N   CYS A  42      40.714  -5.292  12.123  1.00 11.29           N",
        "ATOM    189  CA  CYS A  42      39.736  -5.883  12.911  1.00 10.01           C",
        "ATOM    190  C   CYS A  42      40.339  -6.654  14.087  1.00 22.28           C",
        "ATOM    191  O   CYS A  42      41.181  -7.530  13.859  1.00 13.70           O",
        "ATOM    192  CB  CYS A  42      38.949  -6.825  12.002  1.00  9.67           C",
        "ATOM    193  SG  CYS A  42      37.557  -7.514  12.922  1.00 20.12           S",
    ]
    positions = np.array(
        [
            [40.714, -5.292, 12.123],
            [39.736, -5.883, 12.911],
            [40.339, -6.654, 14.087],
            [41.181, -7.53, 13.859],
            [38.949, -6.825, 12.002],
            [37.557, -7.514, 12.922],
        ],
    )

    res = pdbstructure.Residue("CYS", 42)
    for line in pdb_lines:
        res._add_atom(pdbstructure.Atom(line))

    for i, c in enumerate(res.iter_positions()):
        eq(c, positions[i])


def test_pdbstructure_2():
    atom = pdbstructure.Atom(
        "ATOM   2209  CB  TYR A 299       6.167  22.607  20.046  1.00  8.12           C",
    )
    expected = np.array([6.167, 22.607, 20.046])
    for i, c in enumerate(atom.iter_coordinates()):
        eq(expected[i], c)


def test_pdbstructure_3():
    loc = pdbstructure.Atom.Location(" ", [1, 2, 3], 1.0, 20.0, "XXX")
    expected = [1, 2, 3]
    for i, c in enumerate(loc):
        eq(expected[i], c)


@flaky_pdb_dl
def test_load_pdb_from_url():
    # load pdb from URL
    t1 = load_pdb("https://www.rcsb.org/pdb/files/4ZUO.pdb.gz")
    t2 = load_pdb("https://www.rcsb.org/pdb/files/4ZUO.pdb")
    eq(t1.n_frames, 1)
    eq(t2.n_frames, 1)
    eq(t1.n_atoms, 6200)
    eq(t2.n_atoms, 6200)


@flaky_pdb_dl
def test_load_from_url():
    # load pdb from URL
    t1 = load("https://www.rcsb.org/pdb/files/4ZUO.pdb.gz")
    t2 = load("https://www.rcsb.org/pdb/files/4ZUO.pdb")
    eq(t1.n_frames, 1)
    eq(t2.n_frames, 1)
    eq(t1.n_atoms, 6200)
    eq(t2.n_atoms, 6200)


def test_3nch_conect(get_fn):
    # This has conect entries that use all available digits, good failure case.
    t1 = load_pdb(get_fn("3nch.pdb.gz"))
    top, bonds = t1.top.to_dataframe()
    bonds = {(a, b): 1 for (a, b, _, _) in bonds}
    eq(bonds[19782, 19783], 1)  # Check that last SO4 molecule has right bonds
    eq(bonds[19782, 19784], 1)  # Check that last SO4 molecule has right bonds
    eq(bonds[19782, 19785], 1)  # Check that last SO4 molecule has right bonds
    eq(bonds[19782, 19786], 1)  # Check that last SO4 molecule has right bonds


def test_3nch_serial_resSeq(get_fn):
    # If you use zero-based indexing, this PDB has quite large gaps in residue and atom numbering,
    # so it's a good test case.  See #528
    # Gold standard values obtained via
    # cat 3nch.pdb |grep ATM|tail -n 5
    # HETATM19787  S   SO4 D 804      -4.788  -9.395  22.515  1.00121.87           S
    # HETATM19788  O1  SO4 D 804      -3.815  -9.511  21.425  1.00105.97           O
    # HETATM19789  O2  SO4 D 804      -5.989  -8.733  21.999  1.00116.13           O
    # HETATM19790  O3  SO4 D 804      -5.130 -10.726  23.043  1.00108.74           O
    # HETATM19791  O4  SO4 D 804      -4.210  -8.560  23.575  1.00112.54           O
    t1 = load_pdb(get_fn("3nch.pdb.gz"))
    top, bonds = t1.top.to_dataframe()

    top2 = Topology.from_dataframe(top, bonds)
    eq(t1.top, top2)

    top = top.set_index("serial")  # Index by the actual data in the PDB
    eq(str(top.loc[19791]["name"]), "O4")
    eq(str(top.loc[19787]["name"]), "S")
    eq(str(top.loc[19787]["resName"]), "SO4")
    eq(int(top.loc[19787]["resSeq"]), 804)


def test_1ncw(get_fn):
    load_pdb(get_fn("1ncw.pdb.gz"))


@flaky_pdb_dl
def test_1vii_url_and_gz(get_fn):
    t1 = load_pdb("https://www.rcsb.org/pdb/files/1vii.pdb.gz")
    t2 = load_pdb("https://www.rcsb.org/pdb/files/1vii.pdb")
    t3 = load_pdb(get_fn("1vii.pdb.gz"))
    t4 = load_pdb(get_fn("1vii.pdb"))
    eq(t1.n_frames, 1)
    eq(t1.n_frames, t2.n_frames)
    eq(t1.n_frames, t3.n_frames)
    eq(t1.n_frames, t4.n_frames)

    eq(t1.n_atoms, t2.n_atoms)
    eq(t1.n_atoms, t3.n_atoms)
    eq(t1.n_atoms, t4.n_atoms)


@flaky_pdb_dl
def test_1vii_load_from_mixture(get_fn):
    # load pdb from URL and locally
    t1 = load(["https://www.rcsb.org/pdb/files/1vii.pdb.gz", get_fn("1vii.pdb.gz")])
    t2 = load_pdb(get_fn("1vii.pdb"))
    t3 = load([get_fn("1vii.pdb"), "https://www.rcsb.org/pdb/files/1vii.pdb", get_fn("1vii.pdb")])

    eq(t1.n_frames, 2)
    eq(t2.n_frames, 1)
    eq(t3.n_frames, 3)

    eq(t1.n_atoms, t2.n_atoms)
    eq(t1.n_atoms, t3.n_atoms)


def test_segment_id(get_fn):
    pdb = load_pdb(get_fn("ala_ala_ala.pdb"))
    pdb.save_pdb(temp)
    pdb2 = load_pdb(temp)

    correct_segment_id = "AAL"
    # check that all segment ids are set correctly
    for ridx, r in enumerate(pdb.top.residues):
        assert r.segment_id == correct_segment_id, (
            "residue %i (0-indexed) does not have segment_id set correctly from ala_ala_ala.pdb" % (ridx)
        )

    # check that all segment ids are set correctly after a new pdb file is written
    for ridx, (r1, r2) in enumerate(zip(pdb.top.residues, pdb2.top.residues)):
        assert r1.segment_id == r2.segment_id, (
            f"segment_id of residue {ridx} (0-indexed) in ala_ala_ala.pdb does not "
            "agree with value in after being written out to a new pdb file"
        )


def test_bfactors(get_fn):
    pdb = load_pdb(get_fn("native.pdb"))
    bfactors0 = np.arange(pdb.n_atoms) / 2.0 - 4.0  # (Get some decimals..)

    pdb.save_pdb(temp, bfactors=bfactors0)

    with open(temp) as fh:
        atom_lines = [line for line in fh.readlines() if re.search(r"^ATOM", line)]

    str_bfactors1 = [line[60:66] for line in atom_lines]
    flt_bfactors1 = np.array([float(i) for i in str_bfactors1])

    # check formatting has a space at the beginning and not at the end
    frmt = np.array([(s[0] == " ") and (s[-1] != " ") for s in str_bfactors1])
    assert np.all(frmt)

    # make sure the numbers are actually the same
    eq(bfactors0, flt_bfactors1)


def test_hex(get_fn):
    pdb = load_pdb(get_fn("water_hex.pdb.gz"))
    assert pdb.n_atoms == 100569
    assert pdb.n_residues == 33523
    pdb.save(temp)


def test_dummy_pdb_box_detection(get_fn, recwarn):
    traj = load(get_fn("2koc.pdb"))
    assert len(recwarn) == 1
    w = recwarn.pop(UserWarning)
    assert "Unlikely unit cell" in str(w.message)
    assert traj.unitcell_lengths is None, "Expected dummy box to be deleted"


def test_multichain_load_cycle(get_fn):
    # Issue 1611, make sure that save/load works for more than 1 chain
    pdb = load(get_fn("issue_1611.pdb"))
    bonds = [(bond.atom1.index, bond.atom2.index) for bond in pdb.topology.bonds]
    pdb.save(temp)
    pdb2 = load_pdb(temp)
    bonds2 = [(bond.atom1.index, bond.atom2.index) for bond in pdb2.topology.bonds]
    assert len(bonds) == len(bonds2)


def test_multichain_load_cycle_noter(get_fn):
    # Issue 1611, make sure that save/load works for more than 1 chain
    pdb = load(get_fn("issue_1611.pdb"))
    bonds = [(bond.atom1.index, bond.atom2.index) for bond in pdb.topology.bonds]
    pdb.save(temp, ter=False)
    pdb2 = load_pdb(temp)
    bonds2 = [(bond.atom1.index, bond.atom2.index) for bond in pdb2.topology.bonds]
    assert len(bonds) == len(bonds2)


def test_load_pdb_input_top(get_fn):
    pdb = get_fn("native.pdb")
    p_1 = load_pdb(pdb)

    p_2 = load_pdb(pdb, top=p_1.topology)

    eq(p_1.xyz, p_2.xyz)
    eq(p_1.topology, p_2.topology)


def test_chimera_indexing(get_fn):
    traj = load_pdb(get_fn("chimera_indexing.pdb"))  # this should just not fail

    assert traj.n_atoms == 9
    assert traj.n_residues == 3
    assert traj.topology._atoms[3].serial == 100000
    assert traj.topology._atoms[3].residue.resSeq == 10000
    assert traj.topology.n_bonds == 6
    assert traj.topology._atoms[6].serial == 100003
    assert traj.topology._atoms[6].residue.resSeq == 10001

    test_pos = np.array(
        [
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
        ],
        dtype=np.float32,
    )
    assert np.array_equal(traj._xyz[0], test_pos)


def test_chimera_indexing_skip(get_fn):
    traj = load_pdb(get_fn("chimera_indexing_skip.pdb"))  # this should just not fail

    assert traj.n_atoms == 6
    assert traj.n_residues == 2
    assert traj.topology._atoms[3].serial == 100003
    assert traj.topology._atoms[3].residue.resSeq == 10001
    assert traj.topology.n_bonds == 4

    test_pos = np.array(
        [
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
        ],
        dtype=np.float32,
    )
    assert np.array_equal(traj._xyz[0], test_pos)


def test_vmd_indexing(get_fn):
    traj = load_pdb(get_fn("vmd_indexing.pdb"))  # this should just not fail

    assert traj.n_atoms == 8
    assert traj.n_residues == 6
    assert traj.topology._atoms[1].residue.resSeq == 2710
    assert traj.topology._atoms[3].serial == 100000
    assert traj.topology._atoms[3].residue.resSeq == 10000
    assert traj.topology.n_bonds == 2
    assert traj.topology._atoms[6].serial == 100003
    assert traj.topology._atoms[6].residue.resSeq == 10001
    assert traj.topology._atoms[7].serial == 100465
    assert traj.topology._atoms[7].residue.resSeq == 11000

    test_pos = np.array(
        [
            [10.613, 0.225, 12.764],
            [10.613, 0.225, 12.764],
            [10.613, 0.225, 12.764],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
        ],
        dtype=np.float32,
    )
    assert np.array_equal(traj._xyz[0], test_pos)


def test_overflow_indexing(get_fn):
    traj = load_pdb(get_fn("overflow_indexing.pdb"))  # this should just not fail

    assert traj.n_atoms == 9
    assert traj.n_residues == 3
    assert traj.topology._atoms[3].serial == 100000
    assert traj.topology._atoms[3].residue.resSeq == 10000
    assert traj.topology.n_bonds == 6
    assert traj.topology._atoms[6].serial == 100003
    assert traj.topology._atoms[6].residue.resSeq == 10001

    test_pos = np.array(
        [
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
        ],
        dtype=np.float32,
    )
    assert np.array_equal(traj._xyz[0], test_pos)


def test_blank_indexing(get_fn):
    with pytest.warns(UserWarning):
        traj = load_pdb(get_fn("blank_indexing.pdb"))  # this should just not fail

    assert traj.n_atoms == 8
    assert traj.n_residues == 6
    assert traj.topology._atoms[1].residue.resSeq == 2710
    assert traj.topology._atoms[3].serial == 100000
    assert traj.topology._atoms[3].residue.resSeq == 10000
    assert traj.topology.n_bonds == 2
    assert traj.topology._atoms[6].serial == 100003
    assert traj.topology._atoms[6].residue.resSeq == 10001
    assert traj.topology._atoms[7].serial == 100004
    assert traj.topology._atoms[7].residue.resSeq == 10002

    test_pos = np.array(
        [
            [10.613, 0.225, 12.764],
            [10.613, 0.225, 12.764],
            [10.613, 0.225, 12.764],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
            [10.596, 0.172, 12.686],
            [10.613, 0.225, 12.764],
            [10.629, 0.313, 12.729],
        ],
        dtype=np.float32,
    )
    assert np.array_equal(traj._xyz[0], test_pos)


def test_pdb_charge_read(get_fn):
    """Test that formal charges are read from PDB file"""

    traj = load_pdb(get_fn("1ply_charge.pdb"))

    charges = [atom.formal_charge for atom in traj.topology.atoms if atom.formal_charge is not None]

    # we expect 10 a list of length 10
    assert len(charges) == 10

    # all of the charges should be 1
    assert set(charges) == {1.0}


def test_pdb_charge_write(get_fn):
    """Test that formal charges are written to PDB file"""

    traj = load_pdb(get_fn("1ply_charge.pdb"))

    # write the trajectory to a new file
    traj.save_pdb(temp)

    # read the new file
    traj2 = load_pdb(temp)

    charges = [atom.formal_charge for atom in traj2.topology.atoms if atom.formal_charge is not None]

    # we expect 10 a list of length 10
    assert len(charges) == 10

    # all of the charges should be 1
    assert set(charges) == {1.0}


@pytest.mark.parametrize(
    "bond_orders, ref_orders",
    [
        (False, [None] * 32),
        (True, [None, 2] + [None] * 11 + [2] + [None] * 17 + [1]),
    ],
    ids=[
        "bond_orders=False",
        "bond_orders=True",
    ],
)
def test_ala3_bond_order_read(get_fn, bond_orders, ref_orders):
    """Test that bond orders/types are read properly from CONECT section of PDB file"""

    traj = load_pdb(get_fn("ala_ala_ala.pdb"), bond_orders=bond_orders)

    # Checking bond orders
    bo = [bond.order for bond in traj.top._bonds]
    eq(bo, ref_orders, err_msg="bond orders do not match reference")

    # Checking bond types
    bt = [bond.type for bond in traj.top._bonds]
    bt_from_bo = [float_to_bond_type(order) for order in bo]
    ref_bt = [float_to_bond_type(order) for order in ref_orders]

    # Make sure bond order and bond type match
    eq(bt_from_bo, bt, err_msg="bond types and bond orders don't match")
    # Make sure bond type match reference
    eq(bt, ref_bt, err_msg="bond types do not match reference")


@pytest.mark.parametrize(
    "bond_orders, ref_orders",
    [
        (False, [1] * 12),
        (True, [2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1]),
    ],
    ids=["omit_write_bond_order", "write_bond_order"],
)
def test_bnz_bond_order_conect_write(get_fn, bond_orders, ref_orders):
    """Test that formal charges are written to PDB file"""

    traj = load_pdb(get_fn("bnz.pdb"), bond_orders=True)

    # write the trajectory to a new file
    traj.save_pdb(temp, bond_orders=bond_orders)

    # read the new file, always with bond order info
    traj2 = load_pdb(temp, bond_orders=True)

    bo = [bond.order for bond in traj2.top._bonds]
    eq(bo, ref_orders, err_msg="Inconsistent number of bonds written to pdb")


@pytest.mark.parametrize(
    "bond_orders, ref_orders",
    [
        (False, ["None", None]),
        (True, ["Triple", 3]),
    ],
    ids=["ignore repeated bonds", "read repeated bonds"],
)
def test_bond_order_conect_priority(get_fn, bond_orders, ref_orders):
    """Test that when reading bond orders from pdb, if the bond is listed
    differently in opposite orders, the highest bond order prevails"""

    traj = load_pdb(get_fn("pdb_repeatbond.pdb"), bond_orders=bond_orders)

    bond = [str(traj.top._bonds[0].type), traj.top._bonds[0].order]
    eq(bond, ref_orders, err_msg="Inconsistent bond type, bond order read")
