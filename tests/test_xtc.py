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
import sys

import numpy as np
from mdtraj import io
from mdtraj.formats import XTCTrajectoryFile
from mdtraj.testing import eq
import pytest


@pytest.fixture()
def fn_xtc(get_fn):
    return get_fn('frame0.xtc')


@pytest.fixture()
def pdb(get_fn):
    return get_fn('native.pdb')


strides = (1, 2, 3, 4, 5, 7, 10, 11)


def test_read_chunk1(get_fn, fn_xtc):
    with XTCTrajectoryFile(fn_xtc, 'r', chunk_size_multiplier=0.5) as f:
        xyz, time, step, box = f.read()

    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    assert eq(xyz, iofile['xyz'])
    assert eq(step, iofile['step'])
    assert eq(box, iofile['box'])
    assert eq(time, iofile['time'])


def test_read_stride(get_fn, fn_xtc):
    # read xtc with stride
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    for s in strides:
        with XTCTrajectoryFile(fn_xtc) as f:
            xyz, time, step, box = f.read(stride=s)
        assert eq(xyz, iofile['xyz'][::s])
        assert eq(step, iofile['step'][::s])
        assert eq(box, iofile['box'][::s])
        assert eq(time, iofile['time'][::s])


def test_read_stride_n_frames(get_fn, fn_xtc):
    # read xtc with stride with n_frames
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    for s in strides:
        with XTCTrajectoryFile(fn_xtc) as f:
            xyz, time, step, box = f.read(n_frames=1000, stride=s)
        assert eq(xyz, iofile['xyz'][::s])
        assert eq(step, iofile['step'][::s])
        assert eq(box, iofile['box'][::s])
        assert eq(time, iofile['time'][::s])


def test_read_stride_offsets(get_fn, fn_xtc):
    # read xtc with stride and offsets
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    for s in strides:
        with XTCTrajectoryFile(fn_xtc) as f:
            f.offsets # pre-compute byte offsets between frames
            xyz, time, step, box = f.read(stride=s)
        assert eq(xyz, iofile['xyz'][::s])
        assert eq(step, iofile['step'][::s])
        assert eq(box, iofile['box'][::s])
        assert eq(time, iofile['time'][::s])


def test_read_stride_n_frames_offsets(get_fn, fn_xtc):
    # read xtc with stride with n_frames and offsets
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    for s in strides:
        with XTCTrajectoryFile(fn_xtc) as f:
            f.offsets # pre-compute byte offsets between frames
            xyz, time, step, box = f.read(n_frames=1000, stride=s)
        assert eq(xyz, iofile['xyz'][::s])
        assert eq(step, iofile['step'][::s])
        assert eq(box, iofile['box'][::s])
        assert eq(time, iofile['time'][::s])


def test_read_stride_switching_offsets(get_fn, fn_xtc):
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    with XTCTrajectoryFile(fn_xtc) as f:
        f.offsets  # pre-compute byte offsets between frames
        # read the first 10 frames with stride of 2
        s = 2
        n_frames = 10
        xyz, time, step, box = f.read(n_frames=n_frames, stride=s)
        assert eq(xyz, iofile['xyz'][:n_frames*s:s])
        assert eq(step, iofile['step'][:n_frames*s:s])
        assert eq(box, iofile['box'][:n_frames*s:s])
        assert eq(time, iofile['time'][:n_frames*s:s])
        # now read the rest with stride 3, should start from frame index 8.
        # eg. np.arange(0, n_frames*s + 1, 2)[-1] == 20
        offset = f.tell()
        assert offset == 20
        s = 3
        xyz, time, step, box = f.read(n_frames=None, stride=s)
        assert eq(xyz, iofile['xyz'][offset::s])
        assert eq(step, iofile['step'][offset::s])
        assert eq(box, iofile['box'][offset::s])
        assert eq(time, iofile['time'][offset::s])


def test_read_atomindices_1(get_fn, fn_xtc):
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    with XTCTrajectoryFile(fn_xtc) as f:
        xyz, time, step, box = f.read(atom_indices=[0, 1, 2])
    assert eq(xyz, iofile['xyz'][:, [0, 1, 2]])
    assert eq(step, iofile['step'])
    assert eq(box, iofile['box'])
    assert eq(time, iofile['time'])


def test_read_atomindices_w_stride(get_fn, fn_xtc):
    # test case for bug: https://github.com/mdtraj/mdtraj/issues/1394
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    for stride in strides:
        with XTCTrajectoryFile(fn_xtc) as f:
            xyz, time, step, box = f.read(atom_indices=[0, 1, 2], stride=stride)
        assert eq(xyz, iofile['xyz'][:, [0, 1, 2]][::stride])
        assert eq(step, iofile['step'][::stride])
        assert eq(box, iofile['box'][::stride])
        assert eq(time, iofile['time'][::stride])


def test_read_atomindices_2(get_fn, fn_xtc):
    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    with XTCTrajectoryFile(fn_xtc) as f:
        xyz, time, step, box = f.read(atom_indices=slice(None, None, 2))
    assert eq(xyz, iofile['xyz'][:, ::2])
    assert eq(step, iofile['step'])
    assert eq(box, iofile['box'])
    assert eq(time, iofile['time'])


def test_read_chunk2(get_fn, fn_xtc):
    with XTCTrajectoryFile(fn_xtc, 'r', chunk_size_multiplier=1) as f:
        xyz, time, step, box = f.read()

    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    assert eq(xyz, iofile['xyz'])
    assert eq(step, iofile['step'])
    assert eq(box, iofile['box'])
    assert eq(time, iofile['time'])


def test_read_chunk3(get_fn, fn_xtc):
    with XTCTrajectoryFile(fn_xtc, chunk_size_multiplier=2) as f:
        xyz, time, step, box = f.read(n_frames=100)

    iofile = io.loadh(get_fn('frame0.xtc.h5'), deferred=False)
    assert eq(xyz, iofile['xyz'][:100])
    assert eq(step, iofile['step'][:100])
    assert eq(box, iofile['box'][:100])
    assert eq(time, iofile['time'][:100])


def test_write_0(tmpdir, fn_xtc):
    with XTCTrajectoryFile(fn_xtc) as f:
        xyz = f.read()[0]

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    f = XTCTrajectoryFile(tmpfn, 'w')
    f.write(xyz)
    f.close()

    with XTCTrajectoryFile(tmpfn) as f:
        xyz2, time2, step2, box2 = f.read()
    eq(xyz, xyz2)


def test_write_1(tmpdir):
    xyz = np.asarray(np.around(np.random.randn(100, 10, 3), 3), dtype=np.float32)
    time = np.asarray(np.random.randn(100), dtype=np.float32)
    step = np.arange(100)
    box = np.asarray(np.random.randn(100, 3, 3), dtype=np.float32)

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w') as f:
        f.write(xyz, time=time, step=step, box=box)
    with XTCTrajectoryFile(tmpfn) as f:
        xyz2, time2, step2, box2 = f.read()

    eq(xyz, xyz2)
    eq(time, time2)
    eq(step, step2)
    eq(box, box2)


def test_write_2(tmpdir):
    xyz = np.asarray(np.around(np.random.randn(100, 10, 3), 3), dtype=np.float32)
    time = np.asarray(np.random.randn(100), dtype=np.float32)
    step = np.arange(100)
    box = np.asarray(np.random.randn(100, 3, 3), dtype=np.float32)

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w') as f:
        for i in range(len(xyz)):
            f.write(xyz[i], time=time[i], step=step[i], box=box[i])
    with XTCTrajectoryFile(tmpfn) as f:
        xyz2, time2, step2, box2 = f.read()

    eq(xyz, xyz2)
    eq(time, time2)
    eq(step, step2)
    eq(box, box2)


def test_read_error_0(tmpdir):
    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with pytest.raises(IOError):
        with XTCTrajectoryFile(tmpfn, 'r') as f:
                f.read()


def test_write_error_0(tmpdir):
    xyz = np.asarray(np.random.randn(100, 3, 3), dtype=np.float32)

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w') as f:
        with pytest.raises(ValueError):
            f.read(xyz)


def test_read_error_1():
    with pytest.raises(IOError):
        XTCTrajectoryFile('/tmp/sdfsdfsdf')


def test_read_error_2(get_fn):
    with pytest.raises(IOError):
        XTCTrajectoryFile(get_fn('frame0.dcd')).read()


def test_xtc_write_wierd_0(tmpdir):
    x0 = np.asarray(np.random.randn(100, 3, 3), dtype=np.float32)
    x1 = np.asarray(np.random.randn(100, 9, 3), dtype=np.float32)
    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w') as f:
        f.write(x0)
        with pytest.raises(ValueError):
            f.write(x1)

    xr = XTCTrajectoryFile(tmpfn).read()[0]
    print(xr.shape)


def test_tell(get_fn):
    with XTCTrajectoryFile(get_fn('frame0.xtc')) as f:
        eq(f.tell(), 0)

        f.read(101)
        eq(f.tell(), 101)

        f.read(3)
        eq(f.tell(), 104)


def test_seek(get_fn):
    reference = XTCTrajectoryFile(get_fn('frame0.xtc')).read()[0]
    with XTCTrajectoryFile(get_fn('frame0.xtc')) as f:
        eq(f.tell(), 0)
        eq(f.read(1)[0][0], reference[0])
        eq(f.tell(), 1)

        xyz = f.read(1)[0][0]
        eq(xyz, reference[1])
        eq(f.tell(), 2)

        f.seek(0)
        eq(f.tell(), 0)
        xyz = f.read(1)[0][0]
        eq(f.tell(), 1)
        eq(xyz, reference[0])

        f.seek(5)  # offset array is going to be built
        assert len(f.offsets) == len(reference)
        eq(f.read(1)[0][0], reference[5])
        eq(f.tell(), 6)

        f.seek(-5, 1)
        eq(f.tell(), 1)
        eq(f.read(1)[0][0], reference[1])


def test_seek_natoms9(tmpdir, get_fn):
    # create a xtc file with 9 atoms and seek it.
    with XTCTrajectoryFile(get_fn('frame0.xtc'), 'r') as fh:
        xyz = fh.read()[0][:, :9, :]

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w', force_overwrite=True) as f:
        f.write(xyz)

    with XTCTrajectoryFile(tmpfn, 'r') as f:
        eq(f.read(1)[0].shape, (1, 9, 3))
        eq(f.tell(), 1)
        f.seek(99)
        eq(f.read(1)[0].squeeze(), xyz[99])
        # seek relative
        f.seek(-1, 1)
        eq(f.read(1)[0].squeeze(), xyz[99])

        f.seek(0, 0)
        eq(f.read(1)[0].squeeze(), xyz[0])


def test_seek_out_of_bounds(get_fn):
    with XTCTrajectoryFile(get_fn('frame0.xtc'), 'r') as fh:
        with pytest.raises(IOError):
            fh.seek(10000000)


def test_ragged_1(tmpdir):
    # try first writing no box vectors,, and then adding some
    xyz = np.random.randn(100, 5, 3)
    time = np.random.randn(100)
    box = np.random.randn(100, 3, 3)

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w', force_overwrite=True) as f:
        f.write(xyz)
        with pytest.raises(ValueError):
            f.write(xyz, time, box)


def test_ragged_2(tmpdir):
    # try first writing no box vectors, and then adding some
    xyz = np.random.randn(100, 5, 3)
    time = np.random.randn(100)
    box = np.random.randn(100, 3, 3)

    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w', force_overwrite=True) as f:
        f.write(xyz, time=time, box=box)
        with pytest.raises(ValueError):
            f.write(xyz)


def test_short_traj(tmpdir):
    tmpfn = '{}/traj.xtc'.format(tmpdir)
    with XTCTrajectoryFile(tmpfn, 'w') as f:
        f.write(np.random.uniform(size=(5, 100000, 3)))
    with XTCTrajectoryFile(tmpfn, 'r') as f:
        assert len(f) == 5, len(f)


not_on_win = pytest.mark.skipif(sys.platform.startswith('win'),
                                reason='Can not open file being written again due to file locking.')
@not_on_win
def test_flush(tmpdir):
    tmpfn = '{}/traj.xtc'.format(tmpdir)
    data = np.random.random((5, 100, 3))
    with XTCTrajectoryFile(tmpfn, 'w') as f:
        f.write(data)
        f.flush()
        # note that f is still open, so we can now try to read the contents flushed to disk.
        with XTCTrajectoryFile(tmpfn, 'r') as f2:
            out = f2.read()
        np.testing.assert_allclose(out[0], data, atol=1E-3)
