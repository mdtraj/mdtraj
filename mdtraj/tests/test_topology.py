##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Kyle A. Beauchamp
# Contributors: Robert McGibbon, Matthew Harrigan
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
import tempfile
import mdtraj as md
import numpy as np
from mdtraj.utils.six.moves import cPickle
from mdtraj.utils import import_
from mdtraj.testing import (get_fn, eq, DocStringFormatTester, skipif,
                            assert_raises)

try:
    from simtk.openmm import app
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False


try:
    import pandas as pd
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False


@skipif(not HAVE_OPENMM)
def test_topology_openmm():
    topology = md.load(get_fn('1bpi.pdb')).topology

    # the openmm trajectory doesn't have the distinction
    # between resSeq and index, so if they're out of whack
    # in the openmm version, that cant be preserved
    for residue in topology.residues:
        residue.resSeq = residue.index
    mm = topology.to_openmm()
    assert isinstance(mm, app.Topology)
    topology2 = md.Topology.from_openmm(mm)
    eq(topology, topology2)


@skipif(not HAVE_OPENMM)
def test_topology_openmm_boxes():
    u = import_('simtk.unit')
    traj = md.load(get_fn('1vii_sustiva_water.pdb'))
    mmtop = traj.topology.to_openmm(traj=traj)
    box = mmtop.getUnitCellDimensions() / u.nanometer


@skipif(not HAVE_PANDAS)
def test_topology_pandas():
    topology = md.load(get_fn('native.pdb')).topology
    atoms, bonds = topology.to_dataframe()

    topology2 = md.Topology.from_dataframe(atoms, bonds)
    eq(topology, topology2)

    topology3 = md.Topology.from_dataframe(atoms)  # Make sure you default arguement of None works, see issue #774
    


@skipif(not HAVE_PANDAS)
def test_topology_pandas_TIP4PEW():
    topology = md.load(get_fn('GG-tip4pew.pdb')).topology
    atoms, bonds = topology.to_dataframe()

    topology2 = md.Topology.from_dataframe(atoms, bonds)
    eq(topology, topology2)

def test_topology_numbers():
    topology = md.load(get_fn('1bpi.pdb')).topology
    assert len(list(topology.atoms)) == topology.n_atoms
    assert len(list(topology.residues)) == topology.n_residues
    assert all([topology.atom(i).index == i for i in range(topology.n_atoms)])

@skipif(not HAVE_PANDAS)
def test_topology_unique_elements_bpti():
    traj = md.load(get_fn('bpti.pdb'))
    top, bonds = traj.top.to_dataframe()
    atoms = np.unique(["C", "O", "N", "H", "S"])
    eq(atoms, np.unique(top.element.values))

def test_chain():
    top = md.load(get_fn('bpti.pdb')).topology
    chain = top.chain(0)
    assert chain.n_residues == len(list(chain.residues))

    atoms = list(chain.atoms)
    assert chain.n_atoms == len(atoms)
    for i in range(chain.n_atoms):
        assert atoms[i] == chain.atom(i)

def test_residue():
    top = md.load(get_fn('bpti.pdb')).topology
    residue = top.residue(0)
    assert len(list(residue.atoms)) == residue.n_atoms
    atoms = list(residue.atoms)
    for i in range(residue.n_atoms):
        assert residue.atom(i) == atoms[i]

def test_nonconsective_resSeq():
    t = md.load(get_fn('nonconsecutive_resSeq.pdb'))
    yield lambda : eq(np.array([r.resSeq for r in t.top.residues]), np.array([1, 3, 5]))
    df1 = t.top.to_dataframe()
    df2 = md.Topology.from_dataframe(*df1).to_dataframe()
    yield lambda : eq(df1[0], df2[0])

    # round-trip through a PDB load/save loop
    fd, fname = tempfile.mkstemp(suffix='.pdb')
    os.close(fd)
    t.save(fname)
    t2 = md.load(fname)
    yield lambda : eq(df1[0], t2.top.to_dataframe()[0])
    os.unlink(fname)

def test_pickle():
    # test pickling of topology (bug #391)
    cPickle.loads(cPickle.dumps(md.load(get_fn('bpti.pdb')).topology))

def test_atoms_by_name():
    top = md.load(get_fn('bpti.pdb')).topology

    atoms = list(top.atoms)
    for atom1, atom2 in zip(top.atoms_by_name('CA'), top.chain(0).atoms_by_name('CA')):
        assert atom1 == atom2
        assert atom1 in atoms
        assert atom1.name == 'CA'

    assert len(list(top.atoms_by_name('CA'))) == sum(1 for _ in atoms if _.name == 'CA')
    assert top.residue(15).atom('CA') == [a for a in top.residue(15).atoms if a.name == 'CA'][0]

    assert_raises(KeyError, lambda: top.residue(15).atom('sdfsdsdf'))

def test_select_atom_indices():
    top = md.load(get_fn('native.pdb')).topology

    yield lambda: eq(top.select_atom_indices('alpha'), np.array([8]))
    yield lambda: eq(top.select_atom_indices('minimal'),
                     np.array([4, 5, 6, 8, 10, 14, 15, 16, 18]))

    assert_raises(ValueError, lambda: top.select_atom_indices('sdfsdfsdf'))

@skipif(not HAVE_OPENMM)
def test_top_dataframe_openmm_roundtrip():
    t = md.load(get_fn('2EQQ.pdb'))
    top, bonds = t.top.to_dataframe()
    t.topology = md.Topology.from_dataframe(top, bonds)
    omm_top = t.top.to_openmm()


def test_n_bonds():
    t = md.load(get_fn('2EQQ.pdb'))
    for atom in t.top.atoms:
        if atom.element.symbol == 'H':
            assert atom.n_bonds == 1
        elif atom.element.symbol == 'C':
            assert atom.n_bonds in [3, 4]
        elif atom.element.symbol == 'O':
            assert atom.n_bonds in [1, 2]


def test_load_unknown_topology():
    try:
        md.load(get_fn('frame0.dcd'), top=get_fn('frame0.dcd'))
    except IOError as e:
        # we want to make sure there's a nice error message than includes
        # a list of the supported topology formats.
        assert all(s in str(e) for s in ('.pdb', '.psf', '.prmtop'))
    else:
        assert False  # fail


def test_select_pairs_args():
    traj = md.load(get_fn('tip3p_300K_1ATM.pdb'))
    assert len(traj.top.select_pairs(selection1='name O', selection2='name O')) == 258 * (258 - 1) // 2
    assert (eq(traj.top.select_pairs(selection1="(name O) or (name =~ 'H.*')", selection2="(name O) or (name =~ 'H.*')").sort(),
               traj.top.select_pairs(selection1='all', selection2='all').sort()))
    assert (eq(traj.top.select_pairs(selection1="name O", selection2="name H1").sort(),
               traj.top.select_pairs(selection1="name H1", selection2="name O").sort()))
    assert (eq(traj.top.select_pairs(selection1=range(traj.n_atoms), selection2="(name O) or (name =~ 'H.*')").sort(),
               traj.top.select_pairs(selection1='all', selection2='all').sort()))
