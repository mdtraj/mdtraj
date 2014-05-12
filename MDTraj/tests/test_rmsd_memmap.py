import os
import numpy as np
import mdtraj as md
from mdtraj.testing import get_fn, assert_raises


def test_1():
    # https://github.com/rmcgibbo/mdtraj/issues/438
    try:
        traj = md.load(get_fn('frame0.h5'))
        np.save('temp.npy', traj.xyz)
        traj.xyz = np.load('temp.npy', mmap_mode='r')

        # since traj isn't precentered, this requires centering
        # the coordinates which is done inplace. but that's not possible
        # with mmap_mode = 'r'
        with assert_raises(ValueError):
            md.rmsd(traj, traj, 0)

        # this should work
        traj.xyz = np.load('temp.npy', mmap_mode='c')
        md.rmsd(traj, traj, 0)

    finally:
        os.unlink('temp.npy')


def test_2():
    # https://github.com/rmcgibbo/mdtraj/issues/438
    try:
        traj = md.load(get_fn('frame0.h5'))
        # precenter the coordinates
        traj.center_coordinates()
        traces = traj._rmsd_traces
        np.save('temp.npy', traj.xyz)
        traj.xyz = np.load('temp.npy', mmap_mode='r')
        traj._rmsd_traces = traces

        # this should work, since we don't need to modify the
        # coordinates inplace
        md.rmsd(traj, traj, 0, precentered=True)

    finally:
        os.unlink('temp.npy')