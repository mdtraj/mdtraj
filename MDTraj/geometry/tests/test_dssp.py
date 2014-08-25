import os
import shutil
import itertools
import tempfile
import subprocess
from distutils.spawn import find_executable

import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif

HAVE_DSSP = find_executable('mkdssp')
DSSP_MSG = "This tests required mkdssp to be installed, from http://swift.cmbi.ru.nl/gv/dssp/"
tmpdir = None

def setup():
    global tmpdir
    tmpdir = tempfile.mkdtemp()

def teardown():
    shutil.rmtree(tmpdir)

def call_dssp(traj, frame=0):
    inp = os.path.join(tmpdir, 'temp.pdb')
    out = os.path.join(tmpdir, 'temp.pdb.dssp')
    traj[frame].save(inp)
    cmd = ['mkdssp', '-i', inp, '-o', out]
    subprocess.check_output(' '.join(cmd), shell=True)

    KEY_LINE = '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA'
    with open(out) as f:
        # exaust the first entries
        max(itertools.takewhile(lambda l: not l.startswith(KEY_LINE), f))
        return ''.join([line[16] for line in f])

def assert_(a, b):
    try:
        assert a == b
    except AssertionError:
        print('a', a)
        print('b', b)
        raise


@skipif(not HAVE_DSSP, DSSP_MSG)
def test_1():
    t = md.load(get_fn('1bpi.pdb'))
    assert_(call_dssp(t), md.compute_dssp(t)[0])

@skipif(not HAVE_DSSP, DSSP_MSG)
def test_2():
    t = md.load(get_fn('2EQQ.pdb'))
    for i in range(len(t)):
        yield lambda: assert_(call_dssp(t[i]), md.compute_dssp(t[i])[0])
