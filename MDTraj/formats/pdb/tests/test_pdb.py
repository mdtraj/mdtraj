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
import re, os, tempfile
from mdtraj.formats.pdb import pdbstructure
from mdtraj.formats.pdb.pdbstructure import PdbStructure
from mdtraj.testing import get_fn, eq, raises
from mdtraj import load, load_pdb
from mdtraj.utils import ilen
from mdtraj import Topology

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

def test_pdbstructure_0():
    pdb_lines = [
        "ATOM    188  N   CYS A  42      40.714  -5.292  12.123  1.00 11.29           N  ",
        "ATOM    189  CA  CYS A  42      39.736  -5.883  12.911  1.00 10.01           C  ",
        "ATOM    190  C   CYS A  42      40.339  -6.654  14.087  1.00 22.28           C  ",
        "ATOM    191  O   CYS A  42      41.181  -7.530  13.859  1.00 13.70           O  ",
        "ATOM    192  CB  CYS A  42      38.949  -6.825  12.002  1.00  9.67           C  ",
        "ATOM    193  SG  CYS A  42      37.557  -7.514  12.922  1.00 20.12           S  "
    ]

    res = pdbstructure.Residue("CYS", 42)
    for l in pdb_lines:
        res._add_atom(pdbstructure.Atom(l))
    for i, atom in enumerate(res):
        eq(pdb_lines[i], str(atom))


def test_pdbstructure_1():
    pdb_lines = [
         "ATOM    188  N   CYS A  42      40.714  -5.292  12.123  1.00 11.29           N",
         "ATOM    189  CA  CYS A  42      39.736  -5.883  12.911  1.00 10.01           C",
         "ATOM    190  C   CYS A  42      40.339  -6.654  14.087  1.00 22.28           C",
         "ATOM    191  O   CYS A  42      41.181  -7.530  13.859  1.00 13.70           O",
         "ATOM    192  CB  CYS A  42      38.949  -6.825  12.002  1.00  9.67           C",
         "ATOM    193  SG  CYS A  42      37.557  -7.514  12.922  1.00 20.12           S"
         ]
    positions = np.array([
        [ 40.714,  -5.292,  12.123],
        [ 39.736,  -5.883,  12.911],
        [ 40.339,  -6.654,  14.087],
        [ 41.181,  -7.53,   13.859],
        [ 38.949,  -6.825,  12.002],
        [ 37.557,  -7.514,  12.922]
    ])
    
    res = pdbstructure.Residue("CYS", 42)
    for l in pdb_lines:
        res._add_atom(pdbstructure.Atom(l))

    for i, c in enumerate(res.iter_positions()):
        eq(c, positions[i])

def test_pdbstructure_2():
    atom = pdbstructure.Atom("ATOM   2209  CB  TYR A 299       6.167  22.607  20.046  1.00  8.12           C")
    expected = np.array([6.167, 22.607, 20.046])
    for i, c in enumerate(atom.iter_coordinates()):
        eq(expected[i], c)

def test_pdbstructure_3():
    loc = pdbstructure.Atom.Location(' ', [1,2,3], 1.0, 20.0, "XXX")
    expected = [1, 2, 3]
    for i, c in enumerate(loc):
        eq(expected[i], c)
        
def test_pdb_from_url():
    # load pdb from URL
    t1 = load_pdb('http://www.rcsb.org/pdb/files/4K6Q.pdb.gz')
    t2 = load_pdb('http://www.rcsb.org/pdb/files/4K6Q.pdb')
    eq(t1.n_frames, 1)
    eq(t2.n_frames, 1)
    eq(t1.n_atoms, 2208)
    eq(t2.n_atoms, 2208)

def test_3nch_conect():
    # This has conect entries that use all available digits, good failure case.
    t1 = load_pdb(get_fn('3nch.pdb.gz'))
    top, bonds = t1.top.to_dataframe()
    bonds = dict(((a, b), 1) for (a, b) in bonds)
    eq(bonds[19782, 19783], 1)  # Check that last SO4 molecule has right bonds
    eq(bonds[19782, 19784], 1)  # Check that last SO4 molecule has right bonds
    eq(bonds[19782, 19785], 1)  # Check that last SO4 molecule has right bonds
    eq(bonds[19782, 19786], 1)  # Check that last SO4 molecule has right bonds


def test_3nch_serial_resSeq():
    # If you use zero-based indexing, this PDB has quite large gaps in residue and atom numbering, so it's a good test case.  See #528
    # Gold standard values obtained via
    # cat 3nch.pdb |grep ATM|tail -n 5
    # HETATM19787  S   SO4 D 804      -4.788  -9.395  22.515  1.00121.87           S  
    # HETATM19788  O1  SO4 D 804      -3.815  -9.511  21.425  1.00105.97           O  
    # HETATM19789  O2  SO4 D 804      -5.989  -8.733  21.999  1.00116.13           O  
    # HETATM19790  O3  SO4 D 804      -5.130 -10.726  23.043  1.00108.74           O  
    # HETATM19791  O4  SO4 D 804      -4.210  -8.560  23.575  1.00112.54           O  
    t1 = load_pdb(get_fn('3nch.pdb.gz'))
    top, bonds = t1.top.to_dataframe()
    
    top2 = Topology.from_dataframe(top, bonds)
    eq(t1.top, top2)
    
    top = top.set_index('serial')  # Index by the actual data in the PDB
    eq(str(top.ix[19791]["name"]), "O4")
    eq(str(top.ix[19787]["name"]), "S")
    eq(str(top.ix[19787]["resName"]), "SO4")
    eq(int(top.ix[19787]["resSeq"]), 804)

def test_1ncw():
    t1 = load_pdb(get_fn('1ncw.pdb.gz'))

def test_1vii_url_and_gz():
    t1 = load_pdb('http://www.rcsb.org/pdb/files/1vii.pdb.gz')
    t2 = load_pdb('http://www.rcsb.org/pdb/files/1vii.pdb')
    t3 = load_pdb(get_fn('1vii.pdb.gz'))
    t4 = load_pdb(get_fn('1vii.pdb'))
    eq(t1.n_frames, 1)
    eq(t1.n_frames, t2.n_frames)
    eq(t1.n_frames, t3.n_frames)
    eq(t1.n_frames, t4.n_frames)
    
    eq(t1.n_atoms, t2.n_atoms)
    eq(t1.n_atoms, t3.n_atoms)
    eq(t1.n_atoms, t4.n_atoms)

def test_bfactors():
    pdb = load_pdb(get_fn('native.pdb'))
    bfactors0 = np.arange(pdb.n_atoms) / 2.0 - 4.0 # (Get some decimals..)

    pdb.save_pdb(temp, bfactors=bfactors0)

    with open(temp, 'r') as fh:
        atom_lines = [line for line in fh.readlines() if re.search(r'^ATOM', line)]

    str_bfactors1 = [l[60:66] for l in atom_lines]
    flt_bfactors1 = np.array([float(i) for i in str_bfactors1])

    # check formatting has a space at the beginning and not at the end
    frmt = np.array([(s[0] == ' ') and (s[-1] != ' ') for s in str_bfactors1])
    assert np.all(frmt)
    
    # make sure the numbers are actually the same
    eq(bfactors0, flt_bfactors1)

def test_hex():
   pdb = load_pdb(get_fn('water_hex.pdb.gz'))
   assert pdb.n_atoms == 100569
   assert pdb.n_residues == 33523
   pdb.save(temp)
