##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
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


from __future__ import print_function, division
import numpy as np
import os, tempfile
from mdtraj import topology
from mdtraj.testing import get_fn, eq, raises
from mdtraj import load
from mdtraj.utils import ilen

pdb = get_fn('native.pdb')
fd, temp = tempfile.mkstemp(suffix='.pdb')
os.close(fd)
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_pdbread():
    p = load(pdb)


def test_pdbwrite():
    p = load(pdb)
    p.save(temp)

    r = load(temp)
    eq(p.xyz, r.xyz)


def test_load_multiframe():
    from mdtraj.pdb.pdbstructure import PdbStructure
    with open(get_fn('multiframe.pdb')) as f:
        pdb = PdbStructure(f)
        yield lambda: eq(len(pdb.models), 2)
        yield lambda: eq(len(pdb.models[0].chains), 1)
        yield lambda: eq(len(pdb.models[0].chains[0].residues), 3)
        yield lambda: eq(ilen(pdb.models[0].iter_atoms()), 22)

        yield lambda: eq(len(pdb.models[1].chains), 1)
        yield lambda: eq(len(pdb.models[1].chains[0].residues), 3)
        yield lambda: eq(ilen(pdb.models[1].iter_atoms()), 22)


    t = load(get_fn('multiframe.pdb'))
    yield lambda: eq(t.n_frames, 2)
    yield lambda: eq(t.n_atoms, 22)
    yield lambda: eq(t.xyz[0], t.xyz[1])


def test_4K6Q():
    t = load(get_fn('4K6Q.pdb'))
    eq(t.n_frames, 1)
    eq(t.n_atoms, 2208)

    # this is a random line from the file
    #ATOM   1567  O   LEU A 201      40.239  19.248 -14.530  1.00 33.74           O

    atom = list(t.top.atoms)[1566]
    eq(atom.element.symbol, 'O')
    eq(atom.residue.name, 'LEU')
    eq(atom.index, 1566)
    eq(t.xyz[0, 1566], np.array([40.239, 19.248, -14.530]) / 10)  # converting to NM

    # this is atom 1913 in the PDB
    atom = list(t.top.atoms)[1911]
    eq(atom.name, 'C1')
    eq(atom.residue.name, 'NAG')
    eq([(a1.index, a2.index) for a1, a2 in t.top.bonds if a1.index == 1911 or a2.index == 1911],
       [(1765, 1911), (1911, 1912), (1911, 1922)])

    # that first bond is from a conect record


def test_2EQQ_0():
    # this is an nmr structure with 20 models
    t = load(get_fn('2EQQ.pdb'))
    yield lambda: eq(t.n_frames, 20)
    yield lambda: eq(t.n_atoms, 423)
    yield lambda: eq(ilen(t.top.residues), 28)


def test_1vii_solvated_with_ligand():
    traj = load(get_fn("1vii_sustiva_water.pdb"))
    eq(len(list(traj.top.bonds)), 5124)
    traj.save(temp)
    traj = load(temp)
    eq(len(list(traj.top.bonds)), 5124)

@raises(ValueError)
def test_write_large():
    traj = load(get_fn('native.pdb'))
    traj.xyz.fill(123456789)
    traj.save(temp)


@raises(ValueError)
def test_write_large_2():
    traj = load(get_fn('native.pdb'))
    traj.xyz.fill(-123456789)
    traj.save(temp)
