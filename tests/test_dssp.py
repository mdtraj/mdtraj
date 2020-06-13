import itertools
import os
import subprocess
from distutils.spawn import find_executable

import mdtraj as md
import numpy as np
import pytest

DSSP_MSG = "This test requires mkdssp to be installed, from http://swift.cmbi.ru.nl/gv/dssp/"
needs_dssp = pytest.mark.skipif(not find_executable('mkdssp'), reason=DSSP_MSG)


def call_dssp(dirname, traj, frame=0):
    inp = os.path.join(dirname, 'temp.pdb')
    out = os.path.join(dirname, 'temp.pdb.dssp')
    traj[frame].save(inp)
    cmd = ['mkdssp', '-i', inp, '-o', out]
    subprocess.check_output(' '.join(cmd), shell=True)

    KEY_LINE = '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA'
    with open(out) as f:
        # exaust the first entries
        max(itertools.takewhile(lambda l: not l.startswith(KEY_LINE), f))
        return np.array([line[16] for line in f if line[13] != '!'])


def assert_(a, b):
    try:
        assert np.all(a == b)
    except AssertionError:
        if len(a) != len(b):
            print('Not the same length: %d vs %s' % (len(a), len(b)))
            raise

        for i, (aa, bb) in enumerate(zip(a, b)):
            if aa == bb:
                print("%3d: '%s' '%s'" % (i, aa, bb))
            else:
                print("%3d: '%s' '%s' <-" % (i, aa, bb))
        raise


@needs_dssp
def test_1(get_fn, tmpdir):
    for fn in ['1bpi.pdb', '1vii.pdb', '4ZUO.pdb', '1am7_protein.pdb']:
        t = md.load_pdb(get_fn(fn))
        t = t.atom_slice(t.top.select_atom_indices('minimal'))
        assert_(call_dssp(tmpdir, t), md.compute_dssp(t, simplified=False)[0])


@needs_dssp
def test_2(get_fn, tmpdir):
    t = md.load(get_fn('2EQQ.pdb'))
    for i in range(len(t)):
        assert_(call_dssp(tmpdir, t[i]), md.compute_dssp(t[i], simplified=False)[0])


@needs_dssp
def test_3(tmpdir):
    # 1COY gives a small error, due to a broken chain.
    pdbids = ['1GAI', '6gsv', '2AAC']
    for pdbid in pdbids:
        t = md.load_pdb('http://www.rcsb.org/pdb/files/%s.pdb' % pdbid)
        t = t.atom_slice(t.top.select_atom_indices('minimal'))
        assert_(call_dssp(tmpdir, t), md.compute_dssp(t, simplified=False)[0])


def test_4(get_fn):
    t = md.load_pdb(get_fn('1am7_protein.pdb'))
    a = md.compute_dssp(t, simplified=True)
    b = md.compute_dssp(t, simplified=False)
    assert len(a) == len(b)
    assert len(a[0]) == len(b[0])
    assert list(np.unique(a[0])) == ['C', 'E', 'H']


def test_5(get_fn):
    t = md.load(get_fn('4waters.pdb'))
    a = md.compute_dssp(t, simplified=True)
    b = md.compute_dssp(t, simplified=False)
    ref = np.array([['NA', 'NA', 'NA', 'NA']])

    np.testing.assert_array_equal(a, ref)
    np.testing.assert_array_equal(b, ref)


def test_7(get_fn):
    t = md.load(get_fn('2EQQ.pdb'))
    a = md.compute_dssp(t, simplified=True)
