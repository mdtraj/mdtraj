##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Alexander Yang
# Contributors:
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


import numpy as np
import mdtraj
from mdtraj import io
from mdtraj.testing import eq
import pytest

def test_read(get_fn):
    filename = get_fn('out.gsd')
    traj = mdtraj.load(filename)

    assert len(traj) == 1000
    assert traj.top.n_atoms == 80
    assert traj.top.n_bonds == 70
    assert traj.top.atom(0).name == 'opls_135'
    assert traj.top.atom(1).name == 'opls_140'

def test_read_start(get_fn):
    filename = get_fn('out.gsd')
    traj = mdtraj.load(filename)
    other = mdtraj.load(filename, start=2)

    eq(traj[2].xyz, other[0].xyz)
    eq(traj[2].unitcell_lengths, other[0].unitcell_lengths)
    assert traj.top == other.top

def test_read_frame(get_fn):
    filename = get_fn('out.gsd')
    traj = mdtraj.load(filename)
    other = mdtraj.load(filename, frame=2)

    eq(traj[2].xyz, other[0].xyz)
    eq(traj[2].unitcell_lengths, other[0].unitcell_lengths)
    assert traj.top == other.top

def test_read_stride(get_fn):
    filename = get_fn('out.gsd')
    traj = mdtraj.load(filename)
    other = mdtraj.load(filename, stride=10)

    eq(traj[0].xyz, other[0].xyz)
    eq(traj[0].unitcell_lengths, other[0].unitcell_lengths)

    eq(traj[10].xyz, other[1].xyz)
    eq(traj[10].unitcell_lengths, other[1].unitcell_lengths)

    assert traj.top == other.top

def test_read_variable_top_error(get_fn):
    filename = get_fn('variable_top.gsd')
    with pytest.raises(IOError):
        traj = mdtraj.load(filename)


def test_write(get_fn, tmpdir):
    filename = get_fn('out.gsd')
    traj = mdtraj.load(filename)
    fn = '{}/compare.gsd'.format(tmpdir)
    traj.save(fn)
    other = mdtraj.load(fn)

    assert traj.top == other.top
    eq(other.xyz, traj.xyz)
    eq(other.unitcell_lengths, traj.unitcell_lengths)

def test_write_frame(get_fn, tmpdir):
    filename = get_fn('out.gsd')
    traj = mdtraj.load(filename)
    fn = '{}/compare.gsd'.format(tmpdir)
    traj[2].save(fn)
    other = mdtraj.load(fn)

    assert traj.top == other.top
    eq(traj[2].xyz, other.xyz)
    eq(traj[2].unitcell_lengths, other.unitcell_lengths)

