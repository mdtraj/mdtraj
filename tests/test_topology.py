##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors: Robert McGibbon, Matthew Harrigan, Carlos Xavier Hernandez
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
import pickle
import tempfile

import mdtraj as md
import numpy as np
import pytest
from mdtraj.testing import eq

try:
    from openmm import app
    import openmm.unit as u
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

needs_openmm = pytest.mark.skipif(not HAVE_OPENMM, reason='needs OpenMM')


@needs_openmm
def test_topology_openmm(get_fn):
    topology = md.load(get_fn('1bpi.pdb')).topology
    topology_with_bond_order = md.load(get_fn('imatinib.mol2')).topology

    # the openmm trajectory doesn't have the distinction
    # between resSeq and index, so if they're out of whack
    # in the openmm version, that cant be preserved
    for top in [topology, topology_with_bond_order]:
        for residue in top.residues:
            residue.resSeq = residue.index
        mm = top.to_openmm()
        assert isinstance(mm, app.Topology)
        topology2 = md.Topology.from_openmm(mm)
        eq(top, topology2)


@needs_openmm
def test_topology_openmm_boxes(get_fn):
    traj = md.load(get_fn('1vii_sustiva_water.pdb'))
    mmtop = traj.topology.to_openmm(traj=traj)
    box = mmtop.getUnitCellDimensions() / u.nanometer


def test_topology_pandas(get_fn):
    topology = md.load(get_fn('native.pdb')).topology
    atoms, bonds = topology.to_dataframe()

    topology2 = md.Topology.from_dataframe(atoms, bonds)
    eq(topology, topology2)

    # Make sure default argument of None works, see issue #774
    topology3 = md.Topology.from_dataframe(atoms)


def test_topology_pandas_TIP4PEW(get_fn):
    topology = md.load(get_fn('GG-tip4pew.pdb')).topology
    atoms, bonds = topology.to_dataframe()

    topology2 = md.Topology.from_dataframe(atoms, bonds)
    eq(topology, topology2)


def test_topology_pandas_2residues_same_resSeq(get_fn):
    topology = md.load(get_fn('two_residues_same_resnum.gro')).topology
    atoms, bonds = topology.to_dataframe()

    topology2 = md.Topology.from_dataframe(atoms, bonds)
    eq(topology, topology2)


def test_topology_numbers(get_fn):
    topology = md.load(get_fn('1bpi.pdb')).topology
    assert len(list(topology.atoms)) == topology.n_atoms
    assert len(list(topology.residues)) == topology.n_residues
    assert all([topology.atom(i).index == i for i in range(topology.n_atoms)])


def test_topology_unique_elements_bpti(get_fn):
    traj = md.load(get_fn('bpti.pdb'))
    top, bonds = traj.top.to_dataframe()
    atoms = np.unique(["C", "O", "N", "H", "S"])
    eq(atoms, np.unique(top.element.values))


def test_chain(get_fn):
    top = md.load(get_fn('bpti.pdb')).topology
    chain = top.chain(0)
    assert chain.n_residues == len(list(chain.residues))

    atoms = list(chain.atoms)
    assert chain.n_atoms == len(atoms)
    for i in range(chain.n_atoms):
        assert atoms[i] == chain.atom(i)


def test_residue(get_fn):
    top = md.load(get_fn('bpti.pdb')).topology
    residue = top.residue(0)
    assert len(list(residue.atoms)) == residue.n_atoms
    atoms = list(residue.atoms)
    for i in range(residue.n_atoms):
        assert residue.atom(i) == atoms[i]


def test_segment_id(get_fn):
    top = md.load(get_fn('ala_ala_ala.pdb')).topology
    assert next(top.residues).segment_id == "AAL", "Segment id is not being assigned correctly for ala_ala_ala.psf"
    df = top.to_dataframe()[0]
    assert len(df["segmentID"] == "AAL") == len(
        df), "Segment id is not being assigned correctly to topology data frame ala_ala_ala.psf"


def test_nonconsective_resSeq(get_fn):
    t = md.load(get_fn('nonconsecutive_resSeq.pdb'))
    assert eq(np.array([r.resSeq for r in t.top.residues]), np.array([1, 3, 5]))
    df1 = t.top.to_dataframe()
    df2 = md.Topology.from_dataframe(*df1).to_dataframe()
    assert eq(df1[0], df2[0])

    # round-trip through a PDB load/save loop
    fd, fname = tempfile.mkstemp(suffix='.pdb')
    os.close(fd)
    t.save(fname)
    t2 = md.load(fname)
    assert eq(df1[0], t2.top.to_dataframe()[0])
    os.unlink(fname)


def test_pickle(get_fn):
    # test pickling of topology (bug #391)
    topology_without_bond_order = md.load(get_fn('bpti.pdb')).topology
    topology_with_bond_order = md.load(get_fn('imatinib.mol2')).topology
    for top in [topology_with_bond_order, topology_without_bond_order]:
        loaded_top = pickle.loads(pickle.dumps(top))
        assert loaded_top == top


def test_atoms_by_name(get_fn):
    top = md.load(get_fn('bpti.pdb')).topology

    atoms = list(top.atoms)
    for atom1, atom2 in zip(top.atoms_by_name('CA'), top.chain(0).atoms_by_name('CA')):
        assert atom1 == atom2
        assert atom1 in atoms
        assert atom1.name == 'CA'

    assert len(list(top.atoms_by_name('CA'))) == sum(1 for _ in atoms if _.name == 'CA')
    assert top.residue(15).atom('CA') == [a for a in top.residue(15).atoms if a.name == 'CA'][0]

    with pytest.raises(KeyError):
        top.residue(15).atom('sdfsdf')


def test_select_atom_indices(get_fn):
    top = md.load(get_fn('native.pdb')).topology

    assert eq(top.select_atom_indices('alpha'), np.array([8]))
    assert eq(top.select_atom_indices('minimal'),
              np.array([4, 5, 6, 8, 10, 14, 15, 16, 18]))

    with pytest.raises(ValueError):
        top.select_atom_indices('sdfsdf')


@needs_openmm
def test_top_dataframe_openmm_roundtrip(get_fn):
    t = md.load(get_fn('2EQQ.pdb'))
    top, bonds = t.top.to_dataframe()
    t.topology = md.Topology.from_dataframe(top, bonds)
    omm_top = t.top.to_openmm()


def test_n_bonds(get_fn):
    t = md.load(get_fn('2EQQ.pdb'))
    for atom in t.top.atoms:
        if atom.element.symbol == 'H':
            assert atom.n_bonds == 1
        elif atom.element.symbol == 'C':
            assert atom.n_bonds in [3, 4]
        elif atom.element.symbol == 'O':
            assert atom.n_bonds in [1, 2]


def test_load_unknown_topology(get_fn):
    try:
        md.load(get_fn('frame0.dcd'), top=get_fn('frame0.dcd'))
    except IOError as e:
        # we want to make sure there's a nice error message than includes
        # a list of the supported topology formats.
        assert all(s in str(e) for s in ('.pdb', '.psf', '.prmtop'))
    else:
        assert False  # fail


def test_unique_pairs():
    n = 10
    a = np.arange(n)
    b = np.arange(n, n + n)

    eq(md.Topology._unique_pairs(a, a).sort(), md.Topology._unique_pairs_equal(a).sort())
    eq(md.Topology._unique_pairs(a, b).sort(), md.Topology._unique_pairs_mutually_exclusive(a, b).sort())


def test_select_pairs(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.pdb'))
    select_pairs = traj.top.select_pairs

    assert len(select_pairs(selection1='name O', selection2='name O')) == 258 * (258 - 1) // 2
    assert len(select_pairs(selection1='name H1', selection2='name O')) == 258 * 258

    selections = iter([
        # Equal
        ("(name O) or (name =~ 'H.*')", "(name O) or (name =~ 'H.*')"),
        ('all', 'all'),

        # Exclusive
        ('name O', 'name H1'),
        ('name H1', 'name O'),

        # Overlap
        (range(traj.n_atoms), 'name O'),
        ('all', 'name O')])

    for select1, select2 in selections:
        select3, select4 = next(selections)
        assert eq(select_pairs(selection1=select1, selection2=select2).sort(),
                  select_pairs(selection1=select3, selection2=select4).sort())


def test_to_fasta(get_fn):
    t = md.load(get_fn('2EQQ.pdb'))
    assert t.topology.to_fasta(0) == "ENFSGGCVAGYMRTPDGRCKPTFYQLIT"


def test_subset(get_fn):
    t1 = md.load(get_fn('2EQQ.pdb')).top
    t2 = t1.subset([1, 2, 3])
    assert t2.n_residues == 1

def test_subset_re_index_residues(get_fn):
    t1 = md.load(get_fn('2EQQ.pdb')).top
    t2 = t1.subset(t1.select('resid 0 2'))
    np.testing.assert_array_equal([0, 1], [rr.index for rr in t2.residues])


def test_molecules(get_fn):
    top = md.load(get_fn('4OH9.pdb')).topology
    molecules = top.find_molecules()
    assert sum(len(mol) for mol in molecules) == top.n_atoms
    assert sum(1 for mol in molecules if len(mol) > 1) == 2  # All but two molecules are water


def test_copy_and_hash(get_fn):
    t = md.load(get_fn('traj.h5'))
    t1 = t.topology
    t2 = t.topology.copy()
    assert t1 == t2

    assert hash(tuple(t1._chains)) == hash(tuple(t2._chains))
    assert hash(tuple(t1._atoms)) == hash(tuple(t2._atoms))
    assert hash(tuple(t1._bonds)) == hash(tuple(t2._bonds))
    assert hash(tuple(t1._residues)) == hash(tuple(t2._residues))

    assert hash(t1) == hash(t2)


def test_topology_sliced_residue_indices(get_fn):
    # https://github.com/mdtraj/mdtraj/issues/1585
    full = md.load(get_fn('1bpi.pdb'))
    residues = full.top.select("resid 1 to 10")
    sliced = full.atom_slice(residues)
    idx = [res.index for res in sliced.top.residues][-1]
    assert idx == sliced.top.n_residues-1
    # Now see if this works
    _ = sliced.topology.residue(idx)

def test_topology_join(get_fn):
    top_1 = md.load(get_fn('2EQQ.pdb')).topology
    top_2 = md.load(get_fn('4OH9.pdb')).topology

    out_topology = top_1.join(top_2)

    eq(out_topology.n_atoms, top_1.n_atoms + top_2.n_atoms)
    eq(out_topology.n_residues, top_1.n_residues + top_2.n_residues)
    eq(top_1.atom(0).residue.name, out_topology.atom(0).residue.name)
    eq(top_2.atom(-1).residue.name, out_topology.atom(-1).residue.name)
    eq(top_1.atom(0).element, out_topology.atom(0).element)
    eq(top_2.atom(-1).element, out_topology.atom(-1).element)

def test_topology_join_keep_resSeq(get_fn):
    top_1 = md.load(get_fn('2EQQ.pdb')).topology
    top_2 = md.load(get_fn('4OH9.pdb')).topology

    out_topology_keepId_True = top_1.join(top_2, keep_resSeq=True)
    out_topology_keepId_False = top_1.join(top_2, keep_resSeq=False)

    out_resSeq_keepId_True = [residue.resSeq for residue in out_topology_keepId_True.residues]
    out_resSeq_keepId_False = [residue.resSeq for residue in out_topology_keepId_False.residues]

    expected_resSeq_keepId_True = (
        [residue.resSeq for residue in top_1.residues
            ] + [
                residue.resSeq for residue in top_2.residues])
    expected_resSeq_keepId_False = list(range(1, len(expected_resSeq_keepId_True) + 1))

    eq(out_resSeq_keepId_True, expected_resSeq_keepId_True)
    eq(out_resSeq_keepId_False, expected_resSeq_keepId_False)
