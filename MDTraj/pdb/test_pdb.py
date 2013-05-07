# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

import os, tempfile
from mdtraj.testing import get_fn, eq
from mdtraj import trajectory

pdb = get_fn('native.pdb')
temp = tempfile.mkstemp(suffix='.pdb')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)


def test_pdbread():
    p = trajectory.load(pdb)


def test_pdbwrite():
    p = trajectory.load(pdb)
    p.save(temp)
    
    os.system('cat %s' % temp)
    r = trajectory.load(temp)
    eq(p.xyz, r.xyz)

def test_load_multiframe():
    from mdtraj.pdb.pdbstructure import PdbStructure
    with open(get_fn('multiframe.pdb')) as f:
        pdb = PdbStructure(f)
        yield lambda: eq(len(pdb.models), 2)
        yield lambda: eq(len(pdb.models[0].chains), 1)
        yield lambda: eq(len(pdb.models[0].chains[0].residues), 3)
        yield lambda: eq(len(list(pdb.models[0].iter_atoms())), 22)

        yield lambda: eq(len(pdb.models[1].chains), 1)
        yield lambda: eq(len(pdb.models[1].chains[0].residues), 3)
        yield lambda: eq(len(list(pdb.models[1].iter_atoms())), 22)


    t = trajectory.load(get_fn('multiframe.pdb'))
    yield lambda: eq(t.n_frames, 2)
    yield lambda: eq(t.n_atoms, 22)
    yield lambda: eq(t.xyz[0], t.xyz[1])