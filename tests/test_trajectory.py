##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A Beauchamp, Matthew Harrigan
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
import functools
from pathlib import Path
from mdtraj.testing import eq
import numpy as np
import mdtraj as md
import mdtraj.utils
from mdtraj.utils import six
from mdtraj.core import element
import mdtraj.core.trajectory
import pytest
import mdtraj.formats
from collections import namedtuple

on_win = (sys.platform == 'win32')
on_py3 = (sys.version_info >= (3, 0))

TrajObj = namedtuple('TrajObj', ['fobj', 'fext', 'fn'])
file_objs = [
    (md.formats.NetCDFTrajectoryFile, 'nc'),
    (md.formats.HDF5TrajectoryFile, 'h5'),
    (md.formats.XTCTrajectoryFile, 'xtc'),
    (md.formats.TRRTrajectoryFile, 'trr'),
    (md.formats.DCDTrajectoryFile, 'dcd'),
    (md.formats.MDCRDTrajectoryFile, 'mdcrd'),
    (md.formats.BINPOSTrajectoryFile, 'binpos'),
    (md.formats.DTRTrajectoryFile, 'dtr'),
    (md.formats.XYZTrajectoryFile, 'xyz'),
    (md.formats.XYZTrajectoryFile, 'xyz.gz'),
    (md.formats.LAMMPSTrajectoryFile, 'lammpstrj'),
    (md.formats.TNGTrajectoryFile, 'tng'),
    (md.formats.LH5TrajectoryFile, 'lh5'),
    (md.formats.PDBTrajectoryFile, 'pdb'),
    (md.formats.PDBTrajectoryFile, 'pdb.gz'),
    (md.formats.AmberNetCDFRestartFile, 'ncrst'),
    (md.formats.AmberRestartFile, 'rst7'),
    (md.formats.GroTrajectoryFile, 'gro')
]


@pytest.fixture(params=file_objs, ids=lambda x: x[1])
def ref_traj(request):
    fobj, fext = request.param
    if (on_win and on_py3):
        if fext == 'lh5':
            pytest.skip("lh5 is not supported on Windows + Py3")

    if fext in ['ncrst', 'rst7']:
        pytest.xfail("These reference files don't exist yet")

    return TrajObj(fobj, fext, "frame0.{}".format(fext))


@pytest.fixture(params=file_objs, ids=lambda x: x[1])
def write_traj(request, tmpdir):
    fobj, fext = request.param
    if (on_win and on_py3):
        if fext == 'lh5':
            pytest.skip("lh5 is not supported on Windows + Py3")

    return TrajObj(fobj, fext, "{}/traj.{}".format(tmpdir, fext))


@pytest.fixture
def write_traj_with_box(write_traj):
    if write_traj.fext in ['binpos', 'xyz', 'xyz.gz', 'pdb', 'pdb.gz', 'lh5']:
        pytest.skip("{} does not store box information".format(write_traj.fext))
    else:
        return write_traj


def test_mismatch(get_fn):
    # loading a 22 atoms xtc with a topology that has 2,000 atoms
    # some kind of error should happen
    with pytest.raises(ValueError):
        md.load(get_fn('frame0.xtc'), top=get_fn('4ZUO.pdb'))


def test_box(get_fn):
    t = md.load(get_fn('native.pdb'))
    assert eq(t.unitcell_vectors, None)
    assert eq(t.unitcell_lengths, None)
    assert eq(t.unitcell_angles, None)
    assert eq(t.unitcell_volumes, None)

    t.unitcell_vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).reshape(1, 3, 3)
    assert eq(np.array([1.0, 1.0, 1.0]), t.unitcell_lengths[0])
    assert eq(np.array([90.0, 90.0, 90.0]), t.unitcell_angles[0])
    assert eq(np.array([1.0]), t.unitcell_volumes)


def test_load_pdb_box(get_fn):
    t = md.load(get_fn('native2.pdb'), no_boxchk=True)
    assert eq(t.unitcell_lengths[0], np.array([0.1, 0.2, 0.3]))
    assert eq(t.unitcell_angles[0], np.array([90.0, 90.0, 90.0]))
    assert eq(t.unitcell_vectors[0], np.array([[0.1, 0, 0], [0, 0.2, 0], [0, 0, 0.3]]))


def test_load_pdb_gz(get_fn):
    t = md.load(get_fn('1ncw.pdb.gz'))
    assert eq(t.n_atoms, 3990)


def test_box_load_save(write_traj_with_box, get_fn):
    t = md.load(get_fn('native2.pdb'), no_boxchk=True)
    top = md.load_topology(get_fn('native.pdb'), no_boxchk=True)

    # make sure than through a load/save
    # cycle, the box information is preserved:
    t.save(write_traj_with_box.fn)
    t2 = md.load(write_traj_with_box.fn, top=top)

    assert t.unitcell_vectors is not None
    assert eq(t.xyz, t2.xyz, decimal=3)
    assert eq(t.unitcell_vectors, t2.unitcell_vectors)
    assert eq(t.unitcell_angles, t2.unitcell_angles)
    assert eq(t.unitcell_lengths, t2.unitcell_lengths)


def test_slice(get_fn):
    t = md.load(get_fn('traj.h5'))

    # with copying
    assert eq((t[0:5] + t[5:10]).xyz, t[0:10].xyz)
    assert eq((t[0:5] + t[5:10]).time, t[0:10].time)
    assert eq((t[0:5] + t[5:10]).unitcell_vectors, t[0:10].unitcell_vectors)
    assert eq((t[0:5] + t[5:10]).unitcell_lengths, t[0:10].unitcell_lengths)
    assert eq((t[0:5] + t[5:10]).unitcell_angles, t[0:10].unitcell_angles)

    # without copying (in place)
    assert eq((t.slice(key=range(5), copy=False) + t.slice(key=range(5, 10), copy=False)).xyz,
              t.slice(key=range(10), copy=False).xyz)
    assert eq((t.slice(key=range(5), copy=False) + t.slice(key=range(5, 10), copy=False)).time,
              t.slice(key=range(10), copy=False).time)
    assert eq((t.slice(key=range(5), copy=False) + t.slice(key=range(5, 10), copy=False)).unitcell_vectors,
              t.slice(key=range(10), copy=False).unitcell_vectors)
    assert eq((t.slice(key=range(5), copy=False) + t.slice(key=range(5, 10), copy=False)).unitcell_lengths,
              t.slice(key=range(10), copy=False).unitcell_lengths)
    assert eq((t.slice(key=range(5), copy=False) + t.slice(key=range(5, 10), copy=False)).unitcell_angles,
              t.slice(key=range(10), copy=False).unitcell_angles)


def test_slice2(get_fn):
    t = md.load(get_fn('traj.h5'))
    # with copying
    assert t[0] == t[[0, 1]][0]
    # without copying (in place)
    assert t.slice(key=0, copy=False) == t.slice(key=[0, 1], copy=True)[0]


def has_time_info(fext):
    # Some formats don't save time information
    return fext not in ['dcd', 'binpos', 'pdb', 'pdb.gz', 'xyz', 'xyz.gz', 'lammpstrj', 'lh5', 'mdcrd']


def precision(fext):
    if fext in ['xyz', 'lammpstrj', 'lh5']:
        return 3
    else:
        return 6


def precision2(fext1, fext2):
    return min(precision(fext1), precision(fext2))


def test_read_path(ref_traj, get_fn):
    top = get_fn('native.pdb')
    t = md.load(Path(get_fn(ref_traj.fn)), top=top)


def test_write_path(write_traj, get_fn):
    if write_traj.fext in ('ncrst', 'rst7'):
        pytest.skip("{} can only store 1 frame per file".format(write_traj.fext))
    if write_traj.fext in ('mdcrd'):
        pytest.skip("{} can only store rectilinear boxes".format(write_traj.fext))
    t = md.load(get_fn('traj.h5'))
    if t.unitcell_vectors is None:
        if write_traj.fext in ('dtr', 'lammpstrj'):
            pytest.skip("{} needs to write unitcells".format(write_traj.fext))
    t.save(Path(write_traj.fn))


def test_read_write(ref_traj, write_traj, get_fn):
    if write_traj.fext in ('ncrst', 'rst7'):
        pytest.skip("{} can only store 1 frame per file".format(write_traj.fext))
    if write_traj.fext in ('mdcrd'):
        pytest.skip("{} can only store rectilinear boxes".format(write_traj.fext))

    top = get_fn('native.pdb')
    t = md.load(get_fn(ref_traj.fn), top=top)

    if t.unitcell_vectors is None:
        if write_traj.fext in ('dtr', 'lammpstrj'):
            pytest.skip("{} needs to write unitcells".format(write_traj.fext))

    t.save(write_traj.fn)
    t2 = md.load(write_traj.fn, top=top)
    eq(t.xyz, t2.xyz, decimal=precision2(ref_traj.fext, write_traj.fext))
    if has_time_info(write_traj.fext):
        eq(t.time, t2.time, decimal=3)


def test_load(ref_traj, get_fn):
    nat = md.load(get_fn('native.pdb'))
    num_block = 3
    t0 = md.load(get_fn(ref_traj.fn), top=nat, discard_overlapping_frames=True)
    t1 = md.load(get_fn(ref_traj.fn), top=nat, discard_overlapping_frames=False)
    t2 = md.load([get_fn(ref_traj.fn) for _ in range(num_block)], top=nat, discard_overlapping_frames=False)
    t3 = md.load([get_fn(ref_traj.fn) for _ in range(num_block)], top=nat, discard_overlapping_frames=True)

    # these don't actually overlap, so discard_overlapping_frames should
    # have no effect. the overlap is between the last frame of one and the
    # first frame of the next.
    assert eq(t0.n_frames, t1.n_frames)
    assert eq(t0.n_frames * num_block, t2.n_frames)
    assert eq(t3.n_frames, t2.n_frames)


def test_hdf5(get_fn):
    t = md.load(get_fn('traj.h5'))
    t2 = md.load(get_fn('native.pdb'))
    t3 = md.load(get_fn('traj.h5'), frame=8)

    assert t.topology == t2.topology
    assert eq(t.time, 0.002 * (1 + np.arange(100)))
    assert eq(t.time, 0.002 * (1 + np.arange(100)))
    assert eq(t[8].xyz, t3.xyz)
    assert eq(t[8].time, t3.time)
    assert eq(t[8].unitcell_vectors, t3.unitcell_vectors)


def test_center(get_fn):
    traj = md.load(get_fn('traj.h5'))
    traj.center_coordinates()
    mu = traj.xyz.mean(1)
    mu0 = np.zeros(mu.shape)
    eq(mu0, mu)

    for a in traj.top.atoms:
        # Set all masses equal so we can compare against unweighted result
        a.element = element.hydrogen

    traj.center_coordinates(mass_weighted=True)
    mu2 = traj.xyz.mean(1)
    eq(mu0, mu2)


def test_center_aind(get_fn):
    traj = md.load(get_fn('traj.h5'))
    traj.restrict_atoms(np.arange(0, traj.n_atoms, 2))
    traj.center_coordinates()
    mu = traj.xyz.mean(1)
    mu0 = np.zeros(mu.shape)
    eq(mu0, mu)

    for a in traj.top.atoms:
        # Set all masses equal so we can compare against unweighted result
        a.element = element.hydrogen

    traj.center_coordinates(mass_weighted=True)
    mu2 = traj.xyz.mean(1)
    eq(mu0, mu2)


def test_float_atom_indices_exception(ref_traj, get_fn):
    # Is an informative error message given when you supply floats for atom_indices?
    top = md.load(get_fn('native.pdb')).topology

    try:
        md.load(get_fn(ref_traj.fn), atom_indices=[0.5, 1.3], top=top)
    except ValueError as e:
        if six.PY3:
            assert e.args[0] == 'indices must be of an integer type. float64 is not an integer type'
        else:
            assert e.message == 'indices must be of an integer type. float64 is not an integer type'


def test_restrict_atoms(get_fn):
    traj = md.load(get_fn('traj.h5'))
    time_address = traj.time.ctypes.data

    desired_atom_indices = [0, 1, 2, 5]
    traj.restrict_atoms(desired_atom_indices)
    atom_indices = [a.index for a in traj.top.atoms]
    eq([0, 1, 2, 3], atom_indices)
    eq(traj.xyz.shape[1], 4)
    eq(traj.n_atoms, 4)
    eq(traj.n_residues, 1)
    eq(len(traj.top._bonds), 2)
    eq(traj.n_residues, traj.topology._numResidues)
    eq(traj.n_atoms, traj.topology._numAtoms)
    eq(np.array([a.index for a in traj.topology.atoms]), np.arange(traj.n_atoms))

    # assert that the time field was not copied
    assert traj.time.ctypes.data == time_address


def test_restrict_atoms_not_inplace(get_fn):
    traj = md.load(get_fn('traj.h5'))
    traj_backup = md.load(get_fn('traj.h5'))
    desired_atom_indices = [0, 1, 2, 5]

    sliced = traj.restrict_atoms(desired_atom_indices, inplace=False)

    # make sure the original one was not modified
    eq(traj.xyz, traj_backup.xyz)
    eq(traj.topology, traj_backup.topology)

    eq(list(range(4)), [a.index for a in sliced.top.atoms])
    eq(sliced.xyz.shape[1], 4)
    eq(sliced.n_atoms, 4)
    eq(sliced.n_residues, 1)
    eq(len(sliced.top._bonds), 2)
    eq(sliced.n_residues, sliced.topology._numResidues)
    eq(sliced.n_atoms, sliced.topology._numAtoms)
    eq(np.array([a.index for a in sliced.topology.atoms]), np.arange(sliced.n_atoms))

    # make sure the two don't alias the same memory
    assert traj.time.ctypes.data != sliced.time.ctypes.data
    assert traj.unitcell_angles.ctypes.data != sliced.unitcell_angles.ctypes.data
    assert traj.unitcell_lengths.ctypes.data != sliced.unitcell_lengths.ctypes.data


def test_array_vs_matrix(get_fn):
    top = md.load(get_fn('native.pdb')).topology
    xyz = np.random.randn(1, 22, 3)
    xyz_mat = np.matrix(xyz)
    t1 = md.Trajectory(xyz, top)
    t2 = md.Trajectory(xyz_mat, top)

    eq(t1.xyz, xyz)
    eq(t2.xyz, xyz)


def test_pdb_unitcell_loadsave(tmpdir, get_fn):
    # Make sure that nonstandard unitcell dimensions are saved and loaded
    # correctly with PDB
    tref = md.load(get_fn('native.pdb'))
    tref.unitcell_lengths = 1 + 0.1 * np.random.randn(tref.n_frames, 3)
    tref.unitcell_angles = 90 + 0.0 * np.random.randn(tref.n_frames, 3)
    fn = "{}/x.pdb".format(tmpdir)
    tref.save(fn)

    tnew = md.load(fn)
    eq(tref.unitcell_vectors, tnew.unitcell_vectors, decimal=3)


def test_load_combination(ref_traj, get_fn):
    # Test that the load function's stride and atom_indices work across
    # all trajectory formats

    topology = md.load(get_fn('native.pdb')).topology
    ainds = np.array([a.index for a in topology.atoms if a.element.symbol == 'C'])

    no_kwargs = md.load(get_fn(ref_traj.fn), top=topology)
    strided3 = md.load(get_fn(ref_traj.fn), top=topology, stride=3)
    subset = md.load(get_fn(ref_traj.fn), top=topology, atom_indices=ainds)

    # test 1
    t1 = no_kwargs
    t2 = strided3
    assert eq(t1.xyz[::3], t2.xyz)
    assert eq(t1.time[::3], t2.time)
    if t1.unitcell_vectors is not None:
        assert eq(t1.unitcell_vectors[::3], t2.unitcell_vectors)
    assert eq(t1.topology, t2.topology)

    # test 2
    t1 = no_kwargs
    t2 = subset
    assert eq(t1.xyz[:, ainds, :], t2.xyz)
    assert eq(t1.time, t2.time)
    if t1.unitcell_vectors is not None:
        assert eq(t1.unitcell_vectors, t2.unitcell_vectors)
    assert eq(t1.topology.subset(ainds), t2.topology)


def test_no_topology():
    # We can make trajectories without a topology
    md.Trajectory(xyz=np.random.randn(10, 5, 3), topology=None)


def test_join():
    xyz = np.random.rand(10, 5, 3)
    # overlapping frames
    t1 = md.Trajectory(xyz=xyz[:5], topology=None)
    t2 = md.Trajectory(xyz=xyz[4:], topology=None)

    t3 = t1.join(t2, discard_overlapping_frames=True)
    t4 = t1.join(t2, discard_overlapping_frames=False)
    eq(t3.xyz, xyz)
    eq(len(t4.xyz), 11)
    eq(t4.xyz, np.vstack((xyz[:5], xyz[4:])))


def test_md_join(get_fn):
    fn = get_fn('traj.h5')
    t_ref = md.load(get_fn('frame0.h5'))[:20]
    loaded = md.load(fn, top=t_ref, stride=2)
    iterloaded = md.join(md.iterload(fn, top=t_ref, stride=2, chunk=6))
    eq(loaded.xyz, iterloaded.xyz)
    eq(loaded.time, iterloaded.time)
    eq(loaded.unitcell_angles, iterloaded.unitcell_angles)
    eq(loaded.unitcell_lengths, iterloaded.unitcell_lengths)


def test_stack_1(get_fn):
    t1 = md.load(get_fn('native.pdb'))
    t2 = t1.stack(t1)
    eq(t2.n_atoms, 2 * t1.n_atoms)
    eq(t2.topology._numAtoms, 2 * t1.n_atoms)
    eq(t1.xyz, t2.xyz[:, 0:t1.n_atoms])
    eq(t1.xyz, t2.xyz[:, t1.n_atoms:])


def test_stack_2():
    t1 = md.Trajectory(xyz=np.random.rand(10, 5, 3), topology=None)
    t2 = md.Trajectory(xyz=np.random.rand(10, 6, 3), topology=None)
    t3 = t1.stack(t2)

    eq(t3.xyz[:, :5], t1.xyz)
    eq(t3.xyz[:, 5:], t2.xyz)
    eq(t3.n_atoms, 11)


def test_seek_read_mode(ref_traj, get_fn):
    # Test the seek/tell capacity of the different TrajectoryFile objects in
    # read mode. Basically, we just seek around the files and read different
    # segments, keeping track of our location manually and checking with both
    # tell() and by checking that the right coordinates are actually returned
    fobj = ref_traj.fobj
    fn = ref_traj.fn

    if ref_traj.fobj is md.formats.PDBTrajectoryFile:
        pytest.xfail("PDB Files don't support seeking")
    if ref_traj.fext == 'xyz.gz':
        pytest.xfail("This is broken")
    if ref_traj.fext == 'gro':
        pytest.xfail("This is broken")

    point = 0
    xyz = md.load(get_fn(fn), top=get_fn('native.pdb')).xyz
    length = len(xyz)
    kwargs = {}
    if fobj is md.formats.MDCRDTrajectoryFile:
        kwargs = {'n_atoms': 22}

    with fobj(get_fn(fn), **kwargs) as f:
        for i in range(100):
            r = np.random.rand()
            if r < 0.25:
                offset = np.random.randint(-5, 5)
                if 0 < point + offset < length:
                    point += offset
                    f.seek(offset, 1)
                else:
                    f.seek(0)
                    point = 0
            if r < 0.5:
                offset = np.random.randint(1, 10)
                if point + offset < length:
                    read = f.read(offset)
                    if fobj not in [md.formats.BINPOSTrajectoryFile, md.formats.LH5TrajectoryFile,
                                    md.formats.XYZTrajectoryFile]:
                        read = read[0]
                    readlength = len(read)
                    read = mdtraj.utils.in_units_of(read, f.distance_unit, 'nanometers')
                    eq(xyz[point:point + offset], read)
                    point += readlength
            elif r < 0.75:
                offset = np.random.randint(low=-100, high=0)
                try:
                    f.seek(offset, 2)
                    point = length + offset
                except NotImplementedError:
                    # not all of the *TrajectoryFiles currently support
                    # seeking from the end, so we'll let this pass if they
                    # say that they dont implement this.
                    pass
            else:
                offset = np.random.randint(100)
                f.seek(offset, 0)
                point = offset

            eq(f.tell(), point)


def test_load_frame(ref_traj, get_fn):
    if ref_traj.fobj is md.formats.GroTrajectoryFile:
        pytest.xfail("Gro doesn't implement seek")
    trajectory = md.load(get_fn(ref_traj.fn), top=get_fn('native.pdb'))
    rand = np.random.randint(len(trajectory))
    frame = md.load_frame(get_fn(ref_traj.fn), index=rand, top=get_fn('native.pdb'))

    if ref_traj.fobj is md.formats.DTRTrajectoryFile:
        pytest.xfail("DTR doesn't load a single frame properly")
    eq(trajectory[rand].xyz, frame.xyz)
    eq(trajectory[rand].unitcell_vectors, frame.unitcell_vectors)
    if has_time_info(ref_traj.fext):
        eq(trajectory[rand].time, frame.time)


def test_load_frame_2eqq(get_fn):
    t1 = md.load(get_fn('2EQQ.pdb'))
    r = np.random.randint(len(t1))
    t2 = md.load_frame(get_fn('2EQQ.pdb'), r)
    eq(t1[r].xyz, t2.xyz)


def test_iterload(write_traj, get_fn):
    if write_traj.fext == 'dtr':
        pytest.xfail("This is broken with dtr")
    t_ref = md.load(get_fn('frame0.h5'))[:20]

    if write_traj.fext in ('ncrst', 'rst7'):
        pytest.skip("Only 1 frame per file format")

    t_ref.save(write_traj.fn)

    for stride in [1, 2, 3]:
        loaded = md.load(write_traj.fn, top=t_ref, stride=stride)
        iterloaded = functools.reduce(lambda a, b: a.join(b),
                                      md.iterload(write_traj.fn, top=t_ref, stride=stride, chunk=6))
        eq(loaded.xyz, iterloaded.xyz)
        eq(loaded.time, iterloaded.time)
        eq(loaded.unitcell_angles, iterloaded.unitcell_angles)
        eq(loaded.unitcell_lengths, iterloaded.unitcell_lengths)


def test_iterload_skip(ref_traj, get_fn):
    if ref_traj.fobj is md.formats.PDBTrajectoryFile:
        pytest.xfail("PDB Iterloads an extra frame!!")
    if ref_traj.fobj is md.formats.GroTrajectoryFile:
        pytest.xfail("Not implemented for some reason")
    if ref_traj.fext in ('ncrst', 'rst7'):
        pytest.skip("Only 1 frame per file format")

    top = md.load(get_fn('native.pdb'))
    t_ref = md.load(get_fn(ref_traj.fn), top=top)

    for cs in [0, 1, 11, 100]:
        for skip in [0, 1, 20, 101]:
            t = functools.reduce(lambda a, b: a.join(b),
                                 md.iterload(get_fn(ref_traj.fn), skip=skip, top=top, chunk=cs))
            eq(t_ref.xyz[skip:], t.xyz)
            eq(t_ref.time[skip:], t.time)
            eq(t_ref.topology, t.topology)


def test_iterload_chunk_dcd(get_fn):
    # Makes sure that the actual chunk size yielded by iterload corresponds to the number of
    # frames specified when calling it (for dcd files).
    file = get_fn("alanine-dipeptide-explicit.dcd")
    top = get_fn("alanine-dipeptide-explicit.pdb")

    skip_frames = 3
    frames_chunk = 2

    full = md.load(file, top=top, stride=skip_frames)
    length = len(full)
    

    chunks = []
    for traj_chunk in md.iterload(file, top=top, stride=skip_frames, chunk=frames_chunk):        
        chunks.append(traj_chunk)
    joined = md.join(chunks)
    assert len(full) == len(joined)
    assert eq(full.xyz, joined.xyz)


def test_save_load(write_traj, get_fn):
    # this cycles all the known formats you can save to, and then tries
    # to reload, using just a single-frame file.

    t_ref = md.load(get_fn('native.pdb'))
    t_ref.unitcell_vectors = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]]])

    t_ref.save(write_traj.fn)
    t = md.load(write_traj.fn, top=t_ref.topology)

    eq(t.xyz, t_ref.xyz)
    eq(t.time, t_ref.time)
    if t._have_unitcell:
        eq(t.unitcell_angles, t_ref.unitcell_angles)
        eq(t.unitcell_lengths, t_ref.unitcell_lengths)


def test_force_overwrite(write_traj, get_fn):
    if write_traj.fext == 'dtr':
        pytest.xfail("This is broken with dtr")

    t_ref = md.load(get_fn('native2.pdb'), no_boxchk=True)
    open(write_traj.fn, 'w').close()
    t_ref.save(write_traj.fn, force_overwrite=True)


def test_force_noverwrite(write_traj, get_fn):
    t_ref = md.load(get_fn('native2.pdb'), no_boxchk=True)
    open(write_traj.fn, 'w').close()
    with pytest.raises(IOError):
        t_ref.save(write_traj.fn, force_overwrite=False)


def test_open_and_load(get_fn):
    # These aren't tested in test_length because they don't support length!
    files = ['frame0.mdcrd', '4waters.arc', 'frame0.lammpstrj']

    for file in files:
        if file.endswith('.mdcrd'):
            opened = md.open(get_fn(file), n_atoms=22)
        else:
            opened = md.open(get_fn(file))

        loaded = md.load(get_fn(file), top=get_fn('native.pdb'))


def test_length(get_fn):
    files = ['frame0.nc', 'frame0.h5', 'frame0.xtc', 'frame0.trr',
             'frame0.dcd', '2EQQ.pdb',
             'frame0.binpos', 'frame0.xyz', 'frame0.tng']
    if not (on_win and on_py3):
        files.append('frame0.lh5')

    for file in files:
        opened = md.open(get_fn(file))

        if '.' + file.rsplit('.', 1)[-1] in mdtraj.core.trajectory._TOPOLOGY_EXTS:
            top = file
        else:
            top = 'native.pdb'

        loaded = md.load(get_fn(file), top=get_fn(top))
        assert len(opened) == len(loaded)


def test_unitcell(write_traj, get_fn):
    # make sure that bogus unitcell vecotrs are not saved

    if write_traj.fext in ['rst7', 'ncrst', 'lammpstrj', 'dtr']:
        pytest.xfail('{} seems to need unit vectors'.format(write_traj.fext))

    top = md.load(get_fn('native.pdb')).restrict_atoms(range(5)).topology
    t = md.Trajectory(xyz=np.random.randn(100, 5, 3), topology=top)
    t.save(write_traj.fn)
    assert eq(md.load(write_traj.fn, top=top).unitcell_vectors, None)


def test_chunk0_iterload(get_fn):
    filename = 'frame0.h5'

    trj0 = md.load(get_fn(filename))

    for trj in md.iterload(get_fn(filename), chunk=0):
        pass

    eq(trj0.n_frames, trj.n_frames)


def test_hashing(get_fn):
    frames = [frame for frame in
              md.iterload(get_fn('frame0.xtc'), chunk=1,
                          top=get_fn('native.pdb'))]
    hashes = [hash(frame) for frame in frames]
    # check all frames have a unique hash value
    assert len(hashes) == len(set(hashes))

    # change topology and ensure hash changes too
    top = frames[0].topology
    top.add_bond(top.atom(0), top.atom(1))

    last_frame_hash = hash(frames[0])
    assert last_frame_hash != hashes[-1]

    # test that trajectories without unitcell data can be hashed
    t1 = md.load(get_fn('1bpi.pdb'))
    t2 = md.load(get_fn('1bpi.pdb'))
    assert hash(t1) == hash(t2)
    

def test_smooth(get_fn):
    from scipy.signal import lfilter, lfilter_zi, filtfilt, butter

    pad = 5
    order = 3
    b, a = butter(order, 2.0 / pad)
    zi = lfilter_zi(b, a)

    signal = np.sin(np.arange(100))
    padded = np.r_[signal[pad - 1: 0: -1], signal, signal[-1: -pad: -1]]

    z, _ = lfilter(b, a, padded, zi=zi * padded[0])
    z2, _ = lfilter(b, a, z, zi=zi * z[0])

    output = filtfilt(b, a, padded)
    test = np.loadtxt(get_fn('smooth.txt'))

    eq(output, test)

@pytest.mark.skip(reason="Broken, maybe only on Python 3.11")
def test_image_molecules(get_fn):
    # Load trajectory with periodic box
    t = md.load(get_fn('alanine-dipeptide-explicit.dcd'), top=get_fn('alanine-dipeptide-explicit.pdb'))
    # Image to new trajectory
    t_new = t.image_molecules(inplace=False)
    # Image inplace without making molecules whole
    t.image_molecules(inplace=True, make_whole=False)
    # Image inplace with making molecules whole
    t.image_molecules(inplace=True, make_whole=True)
    # Image with specified anchor molecules
    molecules = t.topology.find_molecules()
    anchor_molecules = molecules[0:3]
    t.image_molecules(inplace=True, anchor_molecules=anchor_molecules)


def test_load_pdb_no_standard_names(get_fn):
    # Minimal test. Standard_names=False will force load_pdb.py
    # to NOT replace any non-standard atom or residue names in the topology
    md.load(get_fn('native2.pdb'), standard_names=False, no_boxchk=True)
    md.load_pdb(get_fn('native2.pdb'), standard_names=False, no_boxchk=True)


def test_load_with_atom_indices(get_fn):
    t1 = md.load(get_fn('frame0.xtc'), top=get_fn('frame0.gro'), atom_indices=[5])
    t2 = md.load(get_fn('frame0.xtc'), top=get_fn('frame0.gro'))
    t2 = t2.atom_slice([5])
    eq(t1.xyz, t2.xyz)
    eq(t1.time, t2.time)


def test_load_with_frame(get_fn):
    t1 = md.load(get_fn('frame0.xtc'), top=get_fn('frame0.pdb'), frame=3)
    t2 = md.load(get_fn('frame0.xtc'), top=get_fn('frame0.pdb'))
    t2 = t2.slice([3])
    eq(t1.xyz, t2.xyz)
    eq(t1.time, t2.time)


def test_add_remove_atoms(get_fn):
    t = md.load(get_fn('aaqaa-wat.pdb'))
    top = t.topology
    old_atoms = list(top.atoms)[:]
    # Add an atom 'MW' at the end of each water molecule
    for r in list(top.residues)[::-1]:
        if r.name != 'HOH': continue
        atoms = list(r.atoms)
        midx = atoms[-1].index+1
        top.insert_atom('MW',None,r,index=midx)
    mwidx = [a.index for a in list(top.atoms) if a.name == 'MW']
    # Check to see whether the 'MW' atoms have the correct index
    assert mwidx == [183+4*i for i in range(83)]
    # Now delete the atoms again
    for r in list(top.residues)[::-1]:
        if r.name != 'HOH': continue
        atoms = list(r.atoms)
        top.delete_atom_by_index(atoms[-1].index)
    roundtrip_atoms = list(top.atoms)[:]
    # Ensure the atoms are the same after a round trip of adding / deleting
    assert old_atoms == roundtrip_atoms
