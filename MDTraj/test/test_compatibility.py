import tempfile, os
import numpy as np

from mdtraj import load
from mdtraj.testing import get_fn, eq
import mdtraj.compatibility


fn = get_fn('legacy_msmbuilder_trj0.lh5')
nat = get_fn('native.pdb')

def test_legacy_hdf1():
    t0 = load(fn, chunk=1)
    t1 = load(fn, chunk=10)
    t2 = load(fn, chunk=100)

    yield lambda: eq(t0.xyz, t1.xyz)
    yield lambda: eq(t0.xyz, t2.xyz)
    yield lambda: t0.topology == load(nat).topology


def test_legacy_hdf2():
    t0 = load(fn, chunk=10, stride=10)
    t1 = load(fn, chunk=20, stride=10)
    t2 = load(fn, chunk=50, stride=10)
    t3 = load(fn, chunk=1, stride=1)

    yield lambda: eq(t0.xyz, t1.xyz)
    yield lambda: eq(t0.xyz, t2.xyz)
    yield lambda: eq(t0.xyz, t3.xyz[::10])
    yield lambda: t0.topology == load(nat).topology
    
def test_legacy_hdf3():
    t0 = load(fn, frame=0)
    t1 = load(fn)

    yield lambda: eq(t0.xyz, t1[0].xyz)
    yield lambda: t0.topology == load(nat).topology