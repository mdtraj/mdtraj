import os
import mdtraj as md
from mdtraj.testing import get_fn, assert_raises
from sklearn.externals.joblib import load, dump


def test_1():
    # https://github.com/rmcgibbo/mdtraj/issues/438
    try:
        traj = md.load(get_fn('frame0.h5'))
        filenames = dump(traj, 'temp')
        traj = load('temp', mmap_mode='r')

        # since traj isn't precentered, this requires centering
        # the coordinates which is done inplace. but that's not possible
        # with mmap_mode = 'r'
        with assert_raises(ValueError):
            md.rmsd(traj, traj, 0)

        # this should work
        traj = load('temp', mmap_mode='c')
        md.rmsd(traj, traj, 0)

    finally:
        for fn in filenames:
            os.unlink(fn)

def test_2():
    # https://github.com/rmcgibbo/mdtraj/issues/438
    try:
        traj = md.load(get_fn('frame0.h5'))
        # precenter the coordinates
        traj.center_coordinates()
        filenames = dump(traj, 'temp')
        traj = load('temp', mmap_mode='r')

        # this should work, since we don't need to modify the
        # coordinates inplace
        md.rmsd(traj, traj, 0, precentered=True)

    finally:
        for fn in filenames:
            os.unlink(fn)
