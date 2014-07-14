from __future__ import print_function
import mdtraj as md
from mdtraj.testing import get_fn, eq, DocStringFormatTester, skipif
from mdtraj.utils import six, fah
import os
import shutil
import tarfile
import tempfile

def test_fah_core17_1():
    filename = get_fn('frame0.xtc')
    tempdir = tempfile.mkdtemp()
    tar_filename = os.path.join(tempdir, "results-000.tar.bz2")
    archive = tarfile.open(tar_filename, mode='w:bz2')    
    
    with tarfile.open(tar_filename, "w:bz2") as tar:
        tar.add(filename, arcname="positions.xtc")
    tar.close()
    
    shutil.copy(tar_filename, os.path.join(tempdir, "results-001.tar.bz2"))

    trj0 = md.load(get_fn("frame0.xtc"), top=get_fn("frame0.h5"))
    output_filename = os.path.join(tempdir, "traj.h5")
    fah.concatenate_core17(tempdir, trj0, output_filename)

    trj = md.load(output_filename)
    eq(trj.n_atoms, trj0.n_atoms)
    eq(trj.n_frames, trj0.n_frames * 2)

    shutil.copy(tar_filename, os.path.join(tempdir, "results-002.tar.bz2"))

    fah.concatenate_core17(tempdir, trj0, output_filename)
    # Should notice the new file and append it to the HDF file.

    trj = md.load(output_filename)
    eq(trj.n_atoms, trj0.n_atoms)
    eq(trj.n_frames, trj0.n_frames * 3)
