import os
import os.path
import numpy as np
import tempfile
import mdtraj as md
import pytest


def test_1(get_fn):
    # https://github.com/mdtraj/mdtraj/issues/438
    try:
        dir = tempfile.mkdtemp()
        fn = os.path.join(dir, 'temp.npy')
        traj = md.load(get_fn('frame0.h5'))
        np.save(fn, traj.xyz)

        traj.xyz = np.load(fn, mmap_mode='r')

        # since traj isn't precentered, this requires centering
        # the coordinates which is done inplace. but that's not possible
        # with mmap_mode = 'r'
        with pytest.raises(ValueError):
            md.rmsd(traj, traj, 0)

        # this should work
        traj.xyz = np.load(fn, mmap_mode='c')
        md.rmsd(traj, traj, 0)

    finally:
        del traj
        os.unlink(fn)
        os.rmdir(dir)


def test_2(get_fn):
    # https://github.com/mdtraj/mdtraj/issues/438
    try:
        dir = tempfile.mkdtemp()
        fn = os.path.join(dir, 'temp.npy')
        traj = md.load(get_fn('frame0.h5'))
        # precenter the coordinates
        traj.center_coordinates()
        traces = traj._rmsd_traces
        np.save(fn, traj.xyz)
        traj.xyz = np.load(fn, mmap_mode='r')
        traj._rmsd_traces = traces

        with pytest.raises(ValueError):
            md.rmsd(traj, traj, 0, precentered=True)

    finally:
        del traj
        os.unlink(fn)
        os.rmdir(dir)
