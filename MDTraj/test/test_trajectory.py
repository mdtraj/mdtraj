# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

import tempfile, os
from mdtraj import binpos, dcd, io
from mdtraj.testing import get_fn, eq, DocStringFormatTester
import numpy as np
from mdtraj.trajectory import load_hdf, load
import mdtraj.trajectory

TestDocstrings = DocStringFormatTester(mdtraj.trajectory, error_on_none=True)

fn = get_fn('frame0.lh5')
nat = get_fn('native.pdb')
temp1 = tempfile.mkstemp(suffix='.xtc')[1]
temp2 = tempfile.mkstemp(suffix='.dcd')[1]
temp3 = tempfile.mkstemp(suffix='.binpos')[1]
temp4 = tempfile.mkstemp(suffix='.trr')[1]
temp5 = tempfile.mkstemp(suffix='.h5')[1]

def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    for e in [temp1, temp2, temp3, temp4]:
        os.unlink(e)

def test_hdf0():
    t0 = load(fn)


def test_box():
    t = load(get_fn('native.pdb'))
    yield lambda: eq(t.unitcell_vectors, None)
    yield lambda: eq(t.unitcell_lengths, None)
    yield lambda: eq(t.unitcell_angles, None)


    t.unitcell_vectors = np.array([[1,0,0], [0,1,0], [0,0,1]]).reshape(1,3,3)
    yield lambda: eq(np.array([1.0, 1.0, 1.0]), t.unitcell_lengths[0])
    yield lambda: eq(np.array([90.0, 90.0, 90.0]), t.unitcell_angles[0])

def test_load_pdb_box():
    t = load(get_fn('native2.pdb'))
    yield lambda: eq(t.unitcell_lengths[0], np.array([0.1, 0.2, 0.3]))
    yield lambda: eq(t.unitcell_angles[0], np.array([90.0, 90.0, 90.0]))
    yield lambda: eq(t.unitcell_vectors[0], np.array([[0.1,0,0], [0,0.2,0], [0,0,0.3]]))


def test_box_load_save():
    t = load(get_fn('native2.pdb'))

    # these three tempfile have extensions (dcd, xtc, trr) that
    # should store the box information. lets make sure than through a load/save
    # cycle, the box information is preserved:
    for temp_fn in [temp1, temp2, temp4, temp5]:
        t.save(temp_fn)
        t2 = load(temp_fn, top=get_fn('native.pdb'))

        assert t.unitcell_vectors != None
        yield lambda: eq(t.xyz, t2.xyz, decimal=3)
        yield lambda: eq(t.unitcell_vectors, t2.unitcell_vectors)
        yield lambda: eq(t.unitcell_angles, t2.unitcell_angles)
        yield lambda: eq(t.unitcell_lengths, t2.unitcell_lengths)



def test_hdf1():
    t0 = load_hdf(fn, top=nat, chunk=1)
    t1 = load_hdf(fn, top=nat, chunk=10)
    t2 = load_hdf(fn, top=nat, chunk=100)

    yield lambda: eq(t0.xyz, t1.xyz)
    yield lambda: eq(t0.xyz, t2.xyz)

def test_hdf2():
    t0 = load_hdf(fn, top=nat, chunk=10, stride=10)
    t1 = load_hdf(fn, top=nat, chunk=20, stride=10)
    t2 = load_hdf(fn, top=nat, chunk=50, stride=10)
    t3 = load_hdf(fn, top=nat, chunk=1, stride=1)

    yield lambda: eq(t0.xyz, t1.xyz)
    yield lambda: eq(t0.xyz, t2.xyz)
    yield lambda: eq(t0.xyz, t3.xyz[::10])

def test_hdf_frame():
    t0 = load_hdf(fn)
    t1 = load_hdf(fn, frame=1)

    yield lambda: eq(t0[1].xyz, t1.xyz)
    yield lambda: eq(t0[1].unitcell_vectors, t1.unitcell_vectors)
    yield lambda: eq(t0[1].unitcell_angles, t1.unitcell_angles)
    yield lambda: eq(t0[1].unitcell_lengths, t1.unitcell_lengths)
    yield lambda: eq(t0[1].time, t1.time)

def test_slice():
    t = load_hdf(fn, top=nat)
    yield lambda: eq((t[0:5] + t[5:10]).xyz, t[0:10].xyz)
    yield lambda: eq((t[0:5] + t[5:10]).time, t[0:10].time)
    yield lambda: eq((t[0:5] + t[5:10]).unitcell_vectors, t[0:10].unitcell_vectors)
    yield lambda: eq((t[0:5] + t[5:10]).unitcell_lengths, t[0:10].unitcell_lengths)
    yield lambda: eq((t[0:5] + t[5:10]).unitcell_angles, t[0:10].unitcell_angles)


def test_slice2():
    t = load_hdf(get_fn('frame1.lh5'))
    yield lambda: t[0] == t[[0,1]][0]

def test_xtc():
    t = mdtraj.trajectory.load(get_fn('frame0.xtc'), top=nat)
    for e in [temp1, temp2, temp3, temp4]:
        def f():
            t.save(e)
            t2 = mdtraj.trajectory.load(e, top=nat)
            eq(t.xyz, t2.xyz, err_msg=e)

            # ony trr and xtc save the time that we read from the original
            # xtc format
            if e.endswith('.trr') or e.endswith('.xtc'):
                eq(t.time, t2.time, err_msg=e)
        yield f

def test_dcd():
    t = mdtraj.trajectory.load(get_fn('frame0.dcd'), top=nat)
    for e in [temp1, temp2, temp3, temp4]:
        def f():
            t.save(e)
            t2 = mdtraj.trajectory.load(e, top=nat)
            eq(t.xyz, t2.xyz, err_msg=e)
            eq(t.time, t2.time, err_msg=e)
        yield f

def test_binpos():
    t = mdtraj.trajectory.load(get_fn('frame0.binpos'), top=nat)
    for e in [temp1, temp2, temp3, temp4]:
        def f():
            t.save(e)
            t2 = mdtraj.trajectory.load(e, top=nat)
            eq(t.xyz, t2.xyz, err_msg=e)
            eq(t.time, t2.time, err_msg=e)
        yield f

def test_load():
    filenames = ["frame0.xtc", "frame0.trr", "frame0.dcd", "frame0.binpos", "frame0.lh5"]
    num_block = 3
    for filename in filenames:
        t0 = mdtraj.trajectory.load(get_fn(filename), top=nat, discard_overlapping_frames=True)
        t1 = mdtraj.trajectory.load(get_fn(filename), top=nat, discard_overlapping_frames=False)
        t2 = mdtraj.trajectory.load([get_fn(filename) for i in xrange(num_block)], top=nat, discard_overlapping_frames=False)
        t3 = mdtraj.trajectory.load([get_fn(filename) for i in xrange(num_block)], top=nat, discard_overlapping_frames=True)
        eq(t0.n_frames, t1.n_frames)
        eq(t0.n_frames * num_block, t2.n_frames)
        eq(t3.n_frames , t0.n_frames * num_block - num_block + 1)
