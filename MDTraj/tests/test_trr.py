from mdtraj import trr, TRRTrajectoryFile, io
import os, tempfile
import numpy as np
from mdtraj.testing import eq, DocStringFormatTester, get_fn
DocStringFormatTester(trr)

temp = tempfile.mkstemp(suffix='.trr')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_1():
    "Write data and read it back"
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


def test_read_stride():
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
         xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
         xyz3, time3, step3, box3, lambd3 = f.read(stride=3)
    yield lambda: eq(xyz[::3], xyz3)
    yield lambda: eq(step[::3], step3)
    yield lambda: eq(box[::3], box3)
    yield lambda: eq(time[::3], time3)

def test_read_stride_2():
    "trr read stride when n_frames is supplied (different path)"
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
         xyz, time, step, box, lambd = f.read()
    with TRRTrajectoryFile(get_fn('frame0.trr')) as f:
         xyz3, time3, step3, box3, lambd3 = f.read(n_frames=1000, stride=3)
    yield lambda: eq(xyz[::3], xyz3)
    yield lambda: eq(step[::3], step3)
    yield lambda: eq(box[::3], box3)
    yield lambda: eq(time[::3], time3)



def test_15():
    "Write data and read it back"
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        f.write(xyz=xyz, time=time, step=step, lambd=lambd)
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read(n_frames=500)

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)


def test_read_atomindices_1():
    with TRRTrajectoryFile(temp) as f:
         xyz, time, step, box, lambd = f.read()
         
    with TRRTrajectoryFile(temp) as f:
         xyz2, time2, step2, box2, lambd2 = f.read(atom_indices=[0,1,2])
    yield lambda: eq(xyz[:, [0,1,2]], xyz2)
    yield lambda: eq(step, step2)
    yield lambda: eq(box, box2)
    yield lambda: eq(lambd, lambd2)
    yield lambda: eq(time, time2)

def test_read_atomindices_2():
    with TRRTrajectoryFile(temp) as f:
         xyz, time, step, box, lambd = f.read()
         
    with TRRTrajectoryFile(temp) as f:
         xyz2, time2, step2, box2, lambd2 = f.read(atom_indices=slice(None, None, 2))
    yield lambda: eq(xyz[:, ::2], xyz2)
    yield lambda: eq(step, step2)
    yield lambda: eq(box, box2)
    yield lambda: eq(lambd, lambd2)
    yield lambda: eq(time, time2)

def test_2():
    """Write data one frame at a time. This checks how the shape is dealt with,
    because each frame is deficient in shape."""
    xyz = np.array(np.random.randn(500,50,3), dtype=np.float32)
    time = np.random.randn(500)
    step = np.arange(500)
    lambd = np.random.randn(500)

    with TRRTrajectoryFile(temp, 'w') as f:
        for i in range(len(xyz)):
            f.write(xyz=xyz[i], time=time[i], step=step[i], lambd=lambd[i])
    with TRRTrajectoryFile(temp) as f:
        xyz2, time2, step2, box2, lambd2 = f.read()

    yield lambda: eq(xyz, xyz2)
    yield lambda: eq(time, time2)
    yield lambda: eq(step, step2)
    yield lambda: eq(lambd, lambd2)

