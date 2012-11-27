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
from mdtraj.trajectory import load_hdf
import mdtraj.trajectory

TestDocstrings = DocStringFormatTester(mdtraj.trajectory, error_on_none=True)

fn = get_fn('frame0.lh5')
nat = get_fn('native.pdb')
temp1 = tempfile.mkstemp(suffix='.xtc')[1]
temp2 = tempfile.mkstemp(suffix='.dcd')[1]
temp3 = tempfile.mkstemp(suffix='.binpos')[1]


def test_hdf1():
    t0 = load_hdf(fn, top=nat, chunk=1)
    t1 = load_hdf(fn, top=nat, chunk=10)
    t2 = load_hdf(fn, top=nat, chunk=100)

    eq(t0.xyz, t1.xyz)
    eq(t0.xyz, t2.xyz)


def test_hdf2():
    t0 = load_hdf(fn, top=nat, chunk=10, stride=10)
    t1 = load_hdf(fn, top=nat, chunk=20, stride=10)
    t2 = load_hdf(fn, top=nat, chunk=50, stride=10)
    t3 = load_hdf(fn, top=nat, chunk=1, stride=1)

    eq(t0.xyz, t1.xyz)
    eq(t0.xyz, t2.xyz)
    eq(t0.xyz, t3.xyz[::10])

def test_slice():
    t = load_hdf(fn, top=nat)
    eq((t[0:5] + t[5:10]).xyz, t[0:10].xyz)
    eq((t[0:5] + t[5:10]).time, t[0:10].time)
    eq((t[0:5] + t[5:10]).box, t[0:10].box)


def test_xtc():
    t = mdtraj.trajectory.load(get_fn('frame0.xtc'), top=nat)
    t.save(temp1)
    t2 = mdtraj.trajectory.load(temp1, top=nat)
    eq(t.xyz, t2.xyz)
    eq(t.time, t2.time)

    t.save(temp2)
    t2 = mdtraj.trajectory.load(temp2, top=nat)
    eq(t.xyz, t2.xyz)

    t.save(temp3)
    t2 = mdtraj.trajectory.load(temp3, top=nat)
    eq(t.xyz, t2.xyz)

def test_dcd():
    t = mdtraj.trajectory.load(get_fn('frame0.dcd'), top=nat)
    t.save(temp1)
    t2 = mdtraj.trajectory.load(temp1, top=nat)
    eq(t.xyz, t2.xyz)
    eq(t.time, t2.time)

    t.save(temp2)
    t2 = mdtraj.trajectory.load(temp2, top=nat)
    eq(t.xyz, t2.xyz)

    t.save(temp3)
    t2 = mdtraj.trajectory.load(temp3, top=nat)
    eq(t.xyz, t2.xyz)

def test_binpos():
    t = mdtraj.trajectory.load(get_fn('frame0.binpos'), top=nat)
    t.save(temp1)
    t2 = mdtraj.trajectory.load(temp1, top=nat)
    eq(t.xyz, t2.xyz)

    t.save(temp2)
    t2 = mdtraj.trajectory.load(temp2, top=nat)
    eq(t.xyz, t2.xyz)

    t.save(temp3)
    t2 = mdtraj.trajectory.load(temp3, top=nat)
    eq(t.xyz, t2.xyz)
