import pytest
import mdtraj as md
from mdtraj.formats import TNGTrajectoryFile
from mdtraj.testing import eq
import os
import tempfile
import numpy as np


def test_load_trajectory(get_fn):
    # Compare a TNG file to the PDB file it was created from.

    pdbtraj = md.load_pdb(get_fn('frame0.pdb'))
    tngtraj = md.load_tng(get_fn('frame0.tng'), top=pdbtraj.topology)
    eq(pdbtraj.n_frames, tngtraj.n_frames)
    eq(pdbtraj.unitcell_vectors, tngtraj.unitcell_vectors)
    eq(pdbtraj.xyz, tngtraj.xyz)


@pytest.mark.skip(reason="Fails on macOS + Python 3.11, need to investigate.")
def test_load_topology(get_fn):
    # Test loading a Topology from a TNG file.

    traj = md.load_tng(get_fn('tng_example.tng'))
    top = traj.topology
    eq(5, top.n_residues)
    eq(10, top.n_bonds)
    bonds = [(a, b) for a, b in top.bonds]
    for res in top.residues:
        eq('HOH', res.name)
        eq(3, res.n_atoms)
        eq(md.element.oxygen, res.atom(0).element)
        eq(md.element.hydrogen, res.atom(1).element)
        eq(md.element.hydrogen, res.atom(2).element)
        assert (res.atom(0), res.atom(1)) in bonds
        assert (res.atom(0), res.atom(2)) in bonds


def test_write(tmpdir):
    # Write a TNG file, then read it back.
    xyz = np.asarray(np.around(np.random.randn(100, 10, 3), 3), dtype=np.float32)
    time = 0.1 * np.arange(0, 100, dtype=np.float32)
    box = np.asarray(np.random.randn(100, 3, 3), dtype=np.float32)

    outfn = "{}/tmp.tng".format(tmpdir)
    with TNGTrajectoryFile(outfn, 'w') as f:
        f.write(xyz, time=time, box=box)
    with TNGTrajectoryFile(outfn) as f:
        xyz2, time2, box2 = f.read()
    eq(xyz, xyz2)
    eq(time, time2)
    eq(box, box2)
