from mdtraj import trr
import os, tempfile
import numpy as np
from mdtraj.testing import get_fn, eq, DocStringFormatTester
DocStringFormatTester(trr)

temp = tempfile.mkstemp(suffix='.trr')[1]

def test1():
    xyz = np.array(np.random.randn(500,10,3), dtype=np.float32)
    
    trr.write(temp, xyz=xyz)
    xyz2 = trr.read(temp, chunk=10)[0]
    
    os.unlink(temp)
    
    eq(xyz, xyz2)
    
def test2():
    xyz = np.array(np.random.randn(500,10,3), dtype=np.float32)

    trr.write(temp, xyz=xyz)
    xyz2 = trr.read(temp, chunk=1000)[0]

    os.unlink(temp)

    eq(xyz, xyz2)