##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
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

from mdtraj import io
from mdtraj.formats import DCDTrajectoryFile, BINPOSTrajectoryFile
from mdtraj.testing import eq

# frame0.binpos was generated from frame0.dcd using
# VMD. The binpos file has one more frame than the dcd file
# because VMD also saved the PDB coordinates
#fn_binpos = get_fn('frame0.binpos')
#fn_dcd = get_fn('frame0.dcd')


def test_read_0(get_fn):
    fn_binpos = get_fn('frame0.binpos')
    fn_dcd = get_fn('frame0.dcd')
    with BINPOSTrajectoryFile(fn_binpos) as f:
        xyz = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2 = f.read()[0]
    xyz3 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')

    assert eq(xyz[1:], xyz2)
    assert eq(xyz, xyz3)


def test_read_stride(get_fn):
    # Read a binpos with stride=3
    fn_binpos = get_fn('frame0.binpos')
    with BINPOSTrajectoryFile(fn_binpos) as f:
        xyz = f.read(stride=3)
    xyz2 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')
    assert eq(xyz, xyz2[::3])


def test_read_stride_2(get_fn):
    # Read a binpos with stride=3 when n_frames is supplied (different code path)
    fn_binpos = get_fn('frame0.binpos')
    with BINPOSTrajectoryFile(fn_binpos) as f:
        xyz = f.read(n_frames=1000, stride=3)
    xyz2 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')
    assert eq(xyz, xyz2[::3])


def test_read_atom_indices(get_fn):
    "Read a binpos with atom_indices as a list"
    fn_binpos = get_fn('frame0.binpos')
    with BINPOSTrajectoryFile(fn_binpos) as f:
        xyz = f.read(atom_indices=[0, 1, 2])
    xyz2 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')
    assert eq(xyz, xyz2[:, [0, 1, 2], :])


def test_read_atom_indices_slice(get_fn):
    "Read a binpos with atom_indices as a slice"
    fn_binpos = get_fn('frame0.binpos')
    with BINPOSTrajectoryFile(fn_binpos) as f:
        xyz = f.read(atom_indices=slice(0, 10, None))
    xyz2 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')
    assert eq(xyz, xyz2[:, 0:10, :])


def test_read_1(get_fn):
    fn_binpos = get_fn('frame0.binpos')
    fn_dcd = get_fn('frame0.dcd')
    with BINPOSTrajectoryFile(fn_binpos, chunk_size_multiplier=0.5) as f:
        xyz = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2 = f.read()[0]
    xyz3 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')

    assert eq(xyz[1:], xyz2)
    assert eq(xyz, xyz3)


def test_read_2(get_fn):
    fn_binpos = get_fn('frame0.binpos')
    fn_dcd = get_fn('frame0.dcd')
    with BINPOSTrajectoryFile(fn_binpos, chunk_size_multiplier=10) as f:
        xyz = f.read()
    with DCDTrajectoryFile(fn_dcd) as f:
        xyz2 = f.read()[0]
    xyz3 = io.loadh(get_fn('frame0.binpos.h5'), 'xyz')

    assert eq(xyz[1:], xyz2)
    assert eq(xyz, xyz3)


def test_write_1(tmpdir):
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    fn = '{}/x.binpos'.format(tmpdir)
    with BINPOSTrajectoryFile(fn, 'w') as f:
        f.write(xyz)

    xyz2 = BINPOSTrajectoryFile(fn).read()
    assert eq(xyz, xyz2)


def test_write_2(tmpdir):
    xyz = np.array(np.random.randn(500, 10, 3), dtype=np.float32)
    fn = '{}/x.binpos'.format(tmpdir)
    with BINPOSTrajectoryFile(fn, 'w') as f:
        f.write(xyz[0:250])
        f.write(xyz[250:])

    xyz2 = BINPOSTrajectoryFile(fn).read()
    assert eq(xyz, xyz2)


def test_tell(get_fn):
    with BINPOSTrajectoryFile(get_fn('frame0.binpos')) as f:
        eq(f.tell(), 0)

        f.read(101)
        eq(f.tell(), 101)

        f.read(3)
        eq(f.tell(), 104)


def test_seek(get_fn):
    reference = BINPOSTrajectoryFile(get_fn('frame0.binpos')).read()
    with BINPOSTrajectoryFile(get_fn('frame0.binpos')) as f:
        xyz = f.read(1)[0]
        eq(xyz, reference[0])
        eq(f.tell(), 1)

        xyz = f.read(1)[0]
        eq(xyz, reference[1])
        eq(f.tell(), 2)

        f.seek(0)
        eq(f.tell(), 0)
        xyz = f.read(1)[0]
        eq(f.tell(), 1)
        eq(xyz, reference[0])

        f.seek(5)
        eq(f.read(1)[0], reference[5])
        eq(f.tell(), 6)

        f.seek(-5, 1)
        eq(f.tell(), 1)
        xyz = f.read(1)[0]
        eq(xyz, reference[1])

        f.seek(0, 2)
        eq(f.tell(), len(reference))
