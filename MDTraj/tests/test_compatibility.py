##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A. Beauchamp
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import os

from mdtraj import load
from mdtraj.testing import get_fn, eq
import mdtraj.compatibility
from mdtraj.utils import ilen
import tempfile


fn = get_fn('legacy_msmbuilder_trj0.lh5')
nat = get_fn('native.pdb')

def test_legacy_hdf1():
    t0 = load(fn, chunk=1)
    t1 = load(fn, chunk=10)
    t2 = load(fn, chunk=100)

    yield lambda: eq(t0.xyz, t1.xyz)
    yield lambda: eq(t0.xyz, t2.xyz)
    yield lambda: t0.topology == load(nat).topology
    yield lambda: eq(ilen(t0.topology.bonds), 14)


def test_legacy_hdf2():
    t0 = load(fn, chunk=10, stride=10)
    t1 = load(fn, chunk=20, stride=10)
    t2 = load(fn, chunk=50, stride=10)
    t3 = load(fn, chunk=1, stride=1)

    yield lambda: eq(t0.xyz, t1.xyz)
    yield lambda: eq(t0.xyz, t2.xyz)
    yield lambda: eq(t0.xyz, t3.xyz[::10])
    yield lambda: t0.topology == load(nat).topology

def test_legacy_hdf3():
    t0 = load(fn, frame=0)
    t1 = load(fn)

    yield lambda: eq(t0.xyz, t1[0].xyz)
    yield lambda: t0.topology == load(nat).topology


def test_legacy_hdf4_save_and_load():
    t0 = load(fn, frame=0)
    t1 = load(fn)

    temp = tempfile.NamedTemporaryFile(suffix=".lh5")
    mdtraj.compatibility.save_legacy_hdf(t0, temp.name)

    t2 = load(temp.name)

    yield lambda: eq(t0.xyz, t2.xyz)
    yield lambda: t0.topology == t2.topology

def test_legacy_hdf5():
    atom_indices = np.arange(10)
    t0 = load(fn, atom_indices=atom_indices)
    eq(t0.n_atoms, 10)
