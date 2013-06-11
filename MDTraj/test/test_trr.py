from mdtraj import trr, TRRTrajectoryFile
import os, tempfile
import numpy as np
from mdtraj.testing import eq, DocStringFormatTester
DocStringFormatTester(trr)

temp = tempfile.mkstemp(suffix='.trr')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test1():
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read()

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)

