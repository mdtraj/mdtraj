import mdtraj as md
from mdtraj.testing import eq

def test_atom_slicing(get_fn):
    t = md.load(get_fn("frame0.xtc"), top=get_fn("frame0.gro"))
    t1 = t.atom_slice([0, 1])
    t2 = t.atom_slice([1, 0])
    eq(t1.xyz, t2.xyz)
    eq(t1.time, t2.time)

def test_atom_slicing_retain_inplace(get_fn):
    t = md.load(get_fn("frame0.xtc"), top=get_fn("frame0.gro"))
    t.atom_slice([0], inplace=True)
    assert t.n_atoms==1
    assert t.top.n_atoms==1


def test_atom_slicing_retain_not_inplace(get_fn):
    t = md.load(get_fn("frame0.xtc"), top=get_fn("frame0.gro"))
    n_atoms = t.n_atoms
    t1 = t.atom_slice([0], inplace=False)
    assert t.n_atoms==n_atoms
    assert t.top.n_atoms==n_atoms
    assert t1.n_atoms==1
    assert t1.top.n_atoms==1

def test_atom_slicing_not_retain_inplace(get_fn):
    t = md.load(get_fn("frame0.xtc"), top=get_fn("frame0.gro"))
    n_atoms = t.n_atoms
    t.atom_slice([0], retain=False, inplace=True)
    assert t.n_atoms==n_atoms-1
    assert t.top.n_atoms==n_atoms-1

def test_atom_slicing_not_retain_not_inplace(get_fn):
    t = md.load(get_fn("frame0.xtc"), top=get_fn("frame0.gro"))
    n_atoms = t.n_atoms
    t1 = t.atom_slice([0], retain=False, inplace=False)
    assert t.n_atoms==n_atoms
    assert t.top.n_atoms==n_atoms
    assert t1.n_atoms==n_atoms-1
    assert t1.top.n_atoms==n_atoms-1