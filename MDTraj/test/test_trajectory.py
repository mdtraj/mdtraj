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

fn = get_fn('frame0.xtc.h5')

def test_hdf1():
    t0 = load_hdf(fn, top=get_fn('native.pdb'), chunk=1)
    t1 = load_hdf(fn, top=get_fn('native.pdb'), chunk=10)
    t2 = load_hdf(fn, top=get_fn('native.pdb'), chunk=100)

    eq(t0.xyz, t1.xyz)
    eq(t0.xyz, t2.xyz)


def test_hdf2():
    t0 = load_hdf(fn, top=get_fn('native.pdb'), chunk=10, stride=10)
    t1 = load_hdf(fn, top=get_fn('native.pdb'), chunk=20, stride=10)
    t2 = load_hdf(fn, top=get_fn('native.pdb'), chunk=50, stride=10)
    t3 = load_hdf(fn, top=get_fn('native.pdb'), chunk=1, stride=1)

    eq(t0.xyz, t1.xyz)
    eq(t0.xyz, t2.xyz)
    eq(t0.xyz, t3.xyz[::10])

def test_slice():
    t = load_hdf(fn, top=get_fn('native.pdb'))
    eq((t[0:5] + t[5:10]).xyz, t[0:10].xyz)
    eq((t[0:5] + t[5:10]).time, t[0:10].time)
    eq((t[0:5] + t[5:10]).box, t[0:10].box)


def test_xtc_1():
    t = mdtraj.trajectory.load(get_fn('frame0.xtc'), top=get_fn('native.pdb'))

def test_dcd_1():
    t = mdtraj.trajectory.load(get_fn('frame0.dcd'), top=get_fn('native.pdb'))

def test_binpos_1():
    t = mdtraj.trajectory.load(get_fn('frame0.binpos'), top=get_fn('native.pdb'))
