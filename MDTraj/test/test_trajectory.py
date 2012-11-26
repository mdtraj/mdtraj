import tempfile, os
from mdtraj import binpos, dcd, io
from mdtraj.testing import get_fn, eq, DocStringFormatTester
import numpy as np
from mdtraj.trajectory import load_hdf
import mdtraj.trajectory

TestDocstrings = DocStringFormatTester(mdtraj.trajectory, error_on_none=True)

fn = get_fn('frame0.xtc.h5')

def test_hdf1():
    t0 = load_hdf(fn, top=get_fn('native.pdb'), chunk=1)
    t1 = load_hdf(fn, top=get_fn('native.pdb'), chunk=10)
    t2 = load_hdf(fn, top=get_fn('native.pdb'), chunk=100)
    
    eq(t0.xyz, t1.xyz)
    eq(t0.xyz, t2.xyz)
    
    
def test_hdf2():
    t0 = load_hdf(fn, top=get_fn('native.pdb'), chunk=10, stride=10)
    t1 = load_hdf(fn, top=get_fn('native.pdb'), chunk=20, stride=10)
    t2 = load_hdf(fn, top=get_fn('native.pdb'), chunk=50, stride=10)
    t3 = load_hdf(fn, top=get_fn('native.pdb'), chunk=1, stride=1)
    
    eq(t0.xyz, t1.xyz)
    eq(t0.xyz, t2.xyz)
    eq(t0.xyz, t3.xyz[::10])

