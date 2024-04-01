import bz2
import gzip

from mdtraj.testing import eq
from mdtraj.utils import open_maybe_zipped


def test_read_gz(tmpdir):
    fn = f"{tmpdir}/read.gz"
    with gzip.GzipFile(fn, "w") as f:
        f.write(b"COOKIE")
    eq(open_maybe_zipped(fn, "r").read(), "COOKIE")


def test_write_gz(tmpdir):
    fn = f"{tmpdir}/write.gz"
    with open_maybe_zipped(fn, "w") as f:
        f.write("COOKIE")
    with gzip.GzipFile(fn, "r") as f:
        eq(f.read().decode("utf-8"), "COOKIE")


def test_read_bz2(tmpdir):
    fn = f"{tmpdir}/read.bz2"
    with bz2.BZ2File(fn, "w") as f:
        f.write(b"COOKIE")
    eq(open_maybe_zipped(fn, "r").read(), "COOKIE")


def test_write_bz2(tmpdir):
    fn = f"{tmpdir}/write.bz2"
    with open_maybe_zipped(fn, "w") as f:
        f.write("COOKIE")
    with bz2.BZ2File(fn, "r") as f:
        eq(f.read().decode("utf-8"), "COOKIE")
