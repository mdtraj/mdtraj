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
    r = trajectory.load(temp)
    eq(p.xyz, r.xyz)
