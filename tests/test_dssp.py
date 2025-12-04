import itertools

import numpy as np
import pytest
from conftest import flaky_pdb_dl

import mdtraj as md


def call_dssp(get_fn, ref_name, frame=0):
    """Read reference output files created with mkdssp v2.2.1"""
    out = get_fn(f"dssp/{ref_name}_{frame}.dssp")

    KEY_LINE = (
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    "
        "O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
    )
    with open(out) as f:
        # exaust the first entries
        max(itertools.takewhile(lambda line: not line.startswith(KEY_LINE), f))
        return np.array([line[16] for line in f if line[13] != "!"], dtype="<U1")


def assert_(a, b):
    try:
        assert np.all(a == b)
    except AssertionError:
        flag = False  # Flag to see if it's a true error or not (e.g. not PPII)

        if len(a) != len(b):
            print("Not the same length: %d vs %s" % (len(a), len(b)))
            raise

        for i, (aa, bb) in enumerate(zip(a, b)):
            if aa == bb:
                print("%3d: '%s' '%s'" % (i, aa, bb))
            elif aa == "P" and bb == " ":
                # Bypassing cases where mkdssp outputs PPII ("P") and mdtraj doesn't.
                print("%3d: '%s' '%s'" % (i, aa, bb))
            else:
                print("%3d: '%s' '%s' <-" % (i, aa, bb))
                flag = True

        if flag:
            raise


@pytest.mark.parametrize("fn", ["1bpi.pdb", "1vii.pdb", "4ZUO.pdb", "1am7_protein.pdb"])
def test_1(get_fn, tmpdir, fn):
    """This test checks dssp assignments for pdb files in tests/data"""
    t = md.load_pdb(get_fn(fn))
    t = t.atom_slice(t.top.select_atom_indices("minimal"))
    assert_(call_dssp(get_fn, fn), md.compute_dssp(t, simplified=False)[0])


@pytest.mark.parametrize("fn", ["2EQQ.pdb"])
def test_2(get_fn, tmpdir, fn):
    """This test checks dssp assignments on different chains for pdb files in tests/data"""
    t = md.load(get_fn(fn))
    for i in range(len(t)):
        assert_(call_dssp(get_fn, fn, i), md.compute_dssp(t[i], simplified=False)[0])


@flaky_pdb_dl
@pytest.mark.parametrize("pdbid", ["1GAI", "6gsv", "2AAC"])
def test_3(get_fn, tmpdir, pdbid):
    """This test checks dssp assignments on pdb files downloaded from rcsb"""
    # 1COY gives a small error, due to a broken chain.
    t = md.load_pdb("https://www.rcsb.org/pdb/files/%s.pdb" % pdbid)
    t = t.atom_slice(t.top.select_atom_indices("minimal"))
    assert_(call_dssp(get_fn, pdbid), md.compute_dssp(t, simplified=False)[0])


def test_4(get_fn):
    t = md.load_pdb(get_fn("1am7_protein.pdb"))
    a = md.compute_dssp(t, simplified=True)
    b = md.compute_dssp(t, simplified=False)
    assert len(a) == len(b)
    assert len(a[0]) == len(b[0])
    assert list(np.unique(a[0])) == ["C", "E", "H"]


def test_5(get_fn):
    t = md.load(get_fn("4waters.pdb"))
    a = md.compute_dssp(t, simplified=True)
    b = md.compute_dssp(t, simplified=False)
    ref = np.array([["NA", "NA", "NA", "NA"]])

    np.testing.assert_array_equal(a, ref)
    np.testing.assert_array_equal(b, ref)


def test_7(get_fn):
    t = md.load(get_fn("2EQQ.pdb"))
    md.compute_dssp(t, simplified=True)
