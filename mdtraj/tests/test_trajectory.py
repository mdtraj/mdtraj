##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A Beauchamp
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


import tempfile, os
import functools
from mdtraj.utils import enter_temp_directory
from mdtraj.testing import get_fn, eq, DocStringFormatTester, assert_raises, SkipTest
import numpy as np
import mdtraj as md
import mdtraj.utils
from mdtraj.utils import six
from mdtraj.utils.six.moves import xrange
from mdtraj.core import element

import mdtraj.core.trajectory
TestDocstrings = DocStringFormatTester(mdtraj.core.trajectory, error_on_none=True)

fn = get_fn('traj.h5')
nat = get_fn('native.pdb')
tmpfns = {}
for suffix, (fd, temp) in {
      'xtc' : tempfile.mkstemp(suffix='.xtc'),
      'dcd' : tempfile.mkstemp(suffix='.dcd'),
      'binpos' : tempfile.mkstemp(suffix='.binpos'),
      'trr' : tempfile.mkstemp(suffix='.trr'),
      'h5' : tempfile.mkstemp(suffix='.h5'),
      'pdb' : tempfile.mkstemp(suffix='.pdb'),
      'pdb.gz' : tempfile.mkstemp(suffix='.pdb.gz'),
      'nc' : tempfile.mkstemp(suffix='.nc'),
      'lh5' : tempfile.mkstemp(suffix='.lh5'),
      'lammpstrj' : tempfile.mkstemp(suffix='.lammpstrj'),
      'xyz' : tempfile.mkstemp(suffix='.xyz')}.items():
    os.close(fd)
    tmpfns[suffix] = temp

def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""

    for e in tmpfns.values():
        os.unlink(e)

def test_mismatch():
    # loading a 22 atoms xtc with a topology that has 2,000 atoms
    # some kind of error should happen!
    assert_raises(ValueError, lambda: md.load(get_fn('frame0.xtc'), top=get_fn('4K6Q.pdb')))


def test_box():
    t = md.load(get_fn('native.pdb'))
    yield lambda: eq(t.unitcell_vectors, None)
    yield lambda: eq(t.unitcell_lengths, None)
    yield lambda: eq(t.unitcell_angles, None)
    yield lambda: eq(t.unitcell_volumes, None)

    t.unitcell_vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).reshape(1, 3, 3)
    yield lambda: eq(np.array([1.0, 1.0, 1.0]), t.unitcell_lengths[0])
    yield lambda: eq(np.array([90.0, 90.0, 90.0]), t.unitcell_angles[0])
    yield lambda: eq(np.array([1.0]), t.unitcell_volumes)


def test_load_pdb_box():
    t = md.load(get_fn('native2.pdb'))
    yield lambda: eq(t.unitcell_lengths[0], np.array([0.1, 0.2, 0.3]))
    yield lambda: eq(t.unitcell_angles[0], np.array([90.0, 90.0, 90.0]))
    yield lambda: eq(t.unitcell_vectors[0], np.array([[0.1, 0, 0], [0, 0.2, 0], [0, 0, 0.3]]))


def test_load_pdb_gz():
    t = md.load(get_fn('1ncw.pdb.gz'))
    yield lambda: eq(t.n_atoms, 3990)


def test_box_load_save():
    t = md.load(get_fn('native2.pdb'))

    # these four tempfile have extensions (dcd, xtc, trr, h5) that
    # should store the box information. lets make sure than through a load/save
    # cycle, the box information is preserved:
    for temp_fn in [tmpfns['xtc'], tmpfns['dcd'], tmpfns['trr'], tmpfns['h5']]:
        t.save(temp_fn)
        if temp_fn.endswith('.h5'):
            t2 = md.load(temp_fn)
        else:
            t2 = md.load(temp_fn, top=get_fn('native.pdb'))

        assert t.unitcell_vectors != None
        yield lambda: eq(t.xyz, t2.xyz, decimal=3)
        yield lambda: eq(t.unitcell_vectors, t2.unitcell_vectors)
        yield lambda: eq(t.unitcell_angles, t2.unitcell_angles)
        yield lambda: eq(t.unitcell_lengths, t2.unitcell_lengths)


def test_slice():
    t = md.load(fn)
    yield lambda: eq((t[0:5] + t[5:10]).xyz, t[0:10].xyz)
    yield lambda: eq((t[0:5] + t[5:10]).time, t[0:10].time)
    yield lambda: eq((t[0:5] + t[5:10]).unitcell_vectors, t[0:10].unitcell_vectors)
    yield lambda: eq((t[0:5] + t[5:10]).unitcell_lengths, t[0:10].unitcell_lengths)
    yield lambda: eq((t[0:5] + t[5:10]).unitcell_angles, t[0:10].unitcell_angles)


def test_slice2():
    t = md.load(get_fn('traj.h5'))
    yield lambda: t[0] == t[[0,1]][0]


def test_xtc():
    t = md.load(get_fn('frame0.xtc'), top=nat)
    for e in [tmpfns['xtc'], tmpfns['dcd'], tmpfns['binpos'], tmpfns['trr'], tmpfns['h5'], tmpfns['pdb'], tmpfns['pdb.gz'], tmpfns['nc']]:
        def f():
            t.save(e)
            t2 = md.load(e, top=nat)
            eq(t.xyz, t2.xyz, err_msg=e)

            # ony trr and xtc save the time that we read from the original
            # xtc format
            if e.endswith('.trr') or e.endswith('.xtc'):
                eq(t.time, t2.time, err_msg=e)
        yield f


def test_dcd():
    t = md.load(get_fn('frame0.dcd'), top=nat)
    for e in [tmpfns['xtc'], tmpfns['dcd'], tmpfns['binpos'], tmpfns['trr'], tmpfns['h5'], tmpfns['pdb'], tmpfns['pdb.gz'], tmpfns['nc']]:
        def f():
            t.save(e)
            t2 = md.load(e, top=nat)
            eq(t.xyz, t2.xyz, err_msg=e)
            eq(t.time, t2.time, err_msg=e)
        yield f

def test_dtr():
    t = md.load(get_fn('ala_dipeptide_trj/clickme.dtr'), top=get_fn('ala_dipeptide.pdb'))
    for e in [tmpfns['xtc'], tmpfns['dcd'], tmpfns['binpos'], tmpfns['trr'], tmpfns['h5'], tmpfns['pdb'], tmpfns['pdb.gz'], tmpfns['nc']]:
        def f():
            t.save(e)
            t2 = md.load(e, top=get_fn('ala_dipeptide.pdb'))

            # change decimal to 3 since the precision is different in different trajectory
            # format
            eq(t.xyz, t2.xyz, decimal=3, err_msg=e)
            #eq(t.time, t2.time, err_msg=e)
        yield f

def test_binpos():
    t = md.load(get_fn('frame0.binpos'), top=nat)
    for e in [tmpfns['xtc'], tmpfns['dcd'], tmpfns['binpos'], tmpfns['trr'], tmpfns['h5'], tmpfns['pdb'], tmpfns['pdb.gz'], tmpfns['nc']]:
        def f():
            t.save(e)
            t2 = md.load(e, top=nat)
            eq(t.xyz, t2.xyz, err_msg=e)
            eq(t.time, t2.time, err_msg=e)
        yield f


def test_load():
    filenames = ["frame0.xtc", "frame0.trr", "frame0.dcd", "frame0.binpos",
                 "traj.h5", 'legacy_msmbuilder_trj0.lh5', 'frame0.nc', six.u('traj.h5'),
                 "frame0.lammpstrj", "frame0.xyz"]
    num_block = 3
    for filename in filenames:
        t0 = md.load(get_fn(filename), top=nat, discard_overlapping_frames=True)
        t1 = md.load(get_fn(filename), top=nat, discard_overlapping_frames=False)
        t2 = md.load([get_fn(filename) for i in xrange(num_block)], top=nat, discard_overlapping_frames=False)
        t3 = md.load([get_fn(filename) for i in xrange(num_block)], top=nat, discard_overlapping_frames=True)

        # these don't actually overlap, so discard_overlapping_frames should have no effect
        # the overlap is between the last frame of one and the first frame of the next.
        yield lambda: eq(t0.n_frames, t1.n_frames)
        yield lambda: eq(t0.n_frames * num_block, t2.n_frames)
        yield lambda: eq(t3.n_frames, t2.n_frames)


def test_hdf5_0():
    t = md.load(get_fn('traj.h5'))
    t2 = md.load(get_fn('native.pdb'))
    t3 = md.load(get_fn('traj.h5'), frame=8)

    print (t.topology, t2.topology)
    assert t.topology == t2.topology
    yield lambda: eq(t.time, 0.002*(1 + np.arange(100)))
    yield lambda: eq(t.time, 0.002*(1 + np.arange(100)))
    yield lambda: eq(t[8].xyz, t3.xyz)
    yield lambda: eq(t[8].time, t3.time)
    yield lambda: eq(t[8].unitcell_vectors, t3.unitcell_vectors)


def test_center():
    traj = md.load(get_fn('traj.h5'))
    traj.center_coordinates()
    mu = traj.xyz.mean(1)
    mu0 = np.zeros(mu.shape)
    eq(mu0, mu)

    for a in traj.top.atoms:
        #a.element.mass = 1.0  # Set all masses equal so we can compare against unweighted result
        a.element = element.hydrogen

    traj.center_coordinates(mass_weighted=True)
    mu2 = traj.xyz.mean(1)
    eq(mu0, mu2)


def test_center_aind():
    traj = md.load(get_fn('traj.h5'))
    traj.restrict_atoms(np.arange(0, traj.n_atoms, 2))
    traj.center_coordinates()
    mu = traj.xyz.mean(1)
    mu0 = np.zeros(mu.shape)
    eq(mu0, mu)

    for a in traj.top.atoms:
        #a.element.mass = 1.0  # Set all masses equal so we can compare against unweighted result
        a.element = element.hydrogen

    traj.center_coordinates(mass_weighted=True)
    mu2 = traj.xyz.mean(1)
    eq(mu0, mu2)


def test_float_atom_indices_exception():
    "Is an informative error message given when you supply floats for atom_indices?"
    top = md.load(get_fn('native.pdb')).topology
    for ext in md._FormatRegistry.loaders.keys():
        try:
            fn = get_fn('frame0' + ext)
        except:
            continue

        try:
            md.load(fn, atom_indices=[0.5, 1.3], top=top)
        except ValueError as e:
            if six.PY3:
                assert e.args[0] == 'indices must be of an integer type. float64 is not an integer type'
            else:
                assert e.message == 'indices must be of an integer type. float64 is not an integer type'
        except Exception as e:
            raise


def test_restrict_atoms():
    traj = md.load(get_fn('traj.h5'))
    time_address = traj.time.ctypes.data

    desired_atom_indices = [0,1,2,5]
    traj.restrict_atoms(desired_atom_indices)
    atom_indices = [a.index for a in traj.top.atoms]
    eq([0,1,2,3], atom_indices)
    eq(traj.xyz.shape[1], 4)
    eq(traj.n_atoms, 4)
    eq(traj.n_residues, 1)
    eq(len(traj.top._bonds), 2)
    eq(traj.n_residues, traj.topology._numResidues)
    eq(traj.n_atoms, traj.topology._numAtoms)
    eq(np.array([a.index for a in traj.topology.atoms]), np.arange(traj.n_atoms))

    # assert that the time field was not copied
    assert traj.time.ctypes.data == time_address


def test_restrict_atoms_not_inplace():
    traj = md.load(get_fn('traj.h5'))
    traj_backup = md.load(get_fn('traj.h5'))
    desired_atom_indices = [0,1,2,5]

    sliced = traj.restrict_atoms(desired_atom_indices, inplace=False)

    # make sure the original one was not modified
    eq(traj.xyz,  traj_backup.xyz)
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


def test_array_vs_matrix():
    top = md.load(get_fn('native.pdb')).topology
    xyz = np.random.randn(1, 22, 3)
    xyz_mat = np.matrix(xyz)
    t1 = md.Trajectory(xyz, top)
    t2 = md.Trajectory(xyz_mat, top)

    eq(t1.xyz, xyz)
    eq(t2.xyz, xyz)

def test_pdb_unitcell_loadsave():
    """Make sure that nonstandard unitcell dimensions are saved and loaded
    correctly with PDB"""
    tref = md.load(get_fn('native.pdb'))
    tref.unitcell_lengths = 1 + 0.1  * np.random.randn(tref.n_frames, 3)
    tref.unitcell_angles = 90 + 0.0  * np.random.randn(tref.n_frames, 3)
    tref.save(tmpfns['pdb'])

    tnew = md.load(tmpfns['pdb'])
    eq(tref.unitcell_vectors, tnew.unitcell_vectors, decimal=3)


def test_load_combination():
    "Test that the load function's stride and atom_indices work across all trajectory formats"

    topology = md.load(get_fn('native.pdb')).topology
    ainds = np.array([a.index for a in topology.atoms if a.element.symbol == 'C'])
    filenames = ['frame0.binpos', 'frame0.dcd', 'frame0.trr', 'frame0.xtc',
                 'frame0.nc', 'frame0.h5', 'frame0.pdb', 'legacy_msmbuilder_trj0.lh5',
                 'frame0.lammpstrj', 'frame0.xyz']

    no_kwargs = [md.load(fn, top=topology) for fn in map(get_fn, filenames)]
    strided3 =  [md.load(fn, top=topology, stride=3) for fn in map(get_fn, filenames)]
    subset =    [md.load(fn, top=topology, atom_indices=ainds) for fn in map(get_fn, filenames)]


    for i, (t1, t2) in enumerate(zip(no_kwargs, strided3)):
        yield lambda: eq(t1.xyz[::3], t2.xyz)
        yield lambda: eq(t1.time[::3], t2.time)
        if t1.unitcell_vectors != None:
            yield lambda: eq(t1.unitcell_vectors[::3], t2.unitcell_vectors)
        yield lambda: eq(t1.topology, t2.topology)

    for i, (t1, t2) in enumerate(zip(no_kwargs, subset)):
        yield lambda: eq(t1.xyz[:, ainds, :], t2.xyz)
        yield lambda: eq(t1.time, t2.time)
        if t1.unitcell_vectors != None:
            yield lambda: eq(t1.unitcell_vectors, t2.unitcell_vectors)
        yield lambda: eq(t1.topology.subset(ainds), t2.topology)

def test_no_topology():
    "We can make trajectories without a topology"
    md.Trajectory(xyz=np.random.randn(10,5,3), topology=None)

def test_join():
    xyz = np.random.rand(10,5,3)
    # overlapping frames
    t1 = md.Trajectory(xyz=xyz[:5], topology=None)
    t2 = md.Trajectory(xyz=xyz[4:], topology=None)

    t3 = t1.join(t2, discard_overlapping_frames=True)
    t4 = t1.join(t2, discard_overlapping_frames=False)
    eq(t3.xyz, xyz)
    eq(len(t4.xyz), 11)
    eq(t4.xyz, np.vstack((xyz[:5], xyz[4:])))

def test_stack_1():
    t1 = md.load(get_fn('native.pdb'))
    t2 = t1.stack(t1)
    eq(t2.n_atoms, 2*t1.n_atoms)
    eq(t2.topology._numAtoms, 2*t1.n_atoms)
    eq(t1.xyz, t2.xyz[:, 0:t1.n_atoms])
    eq(t1.xyz, t2.xyz[:, t1.n_atoms:])


def test_stack_2():
    t1 = md.Trajectory(xyz=np.random.rand(10,5,3), topology=None)
    t2 = md.Trajectory(xyz=np.random.rand(10,6,3), topology=None)
    t3 = t1.stack(t2)

    eq(t3.xyz[:, :5], t1.xyz)
    eq(t3.xyz[:, 5:], t2.xyz)
    eq(t3.n_atoms, 11)


def test_seek_read_mode():
    # Test the seek/tell capacity of the different TrajectoryFile objects in
    # read mode. Basically, we just seek around the files and read different
    # segments, keeping track of our location manually and checking with both
    # tell() and by checking that the right coordinates are actually returned
    files = [(md.formats.NetCDFTrajectoryFile, 'frame0.nc'),
             (md.formats.HDF5TrajectoryFile, 'frame0.h5'),
             (md.formats.XTCTrajectoryFile, 'frame0.xtc'),
             (md.formats.TRRTrajectoryFile, 'frame0.trr'),
             (md.formats.DCDTrajectoryFile, 'frame0.dcd'),
             (md.formats.MDCRDTrajectoryFile, 'frame0.mdcrd'),
             (md.formats.BINPOSTrajectoryFile, 'frame0.binpos'),
             (md.formats.LH5TrajectoryFile, 'legacy_msmbuilder_trj0.lh5'),
             (md.formats.DTRTrajectoryFile, 'frame0.dtr/clickme.dtr'),
             (md.formats.XYZTrajectoryFile, 'frame0.xyz'),
             (md.formats.LAMMPSTrajectoryFile, 'frame0.lammpstrj'),]

    for a, b in files:
        point = 0
        xyz = md.load(get_fn(b), top=get_fn('native.pdb')).xyz
        length = len(xyz)
        kwargs = {}
        if a is md.formats.MDCRDTrajectoryFile:
            kwargs = {'n_atoms': 22}

        with a(get_fn(b), **kwargs) as f:
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
                        if a not in [md.formats.BINPOSTrajectoryFile, md.formats.LH5TrajectoryFile,
                                     md.formats.XYZTrajectoryFile]:
                            read = read[0]
                        readlength = len(read)
                        read = mdtraj.utils.in_units_of(read, f.distance_unit, 'nanometers')
                        eq(xyz[point:point+offset], read)
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


def test_load_frame():
    files = ['frame0.nc', 'frame0.h5', 'frame0.xtc', 'frame0.trr',
             'frame0.dcd', 'frame0.mdcrd', 'frame0.binpos',
             'legacy_msmbuilder_trj0.lh5', 'frame0.xyz', 'frame0.lammpstrj']
    trajectories = [md.load(get_fn(f), top=get_fn('native.pdb')) for f in files]
    rand = [np.random.randint(len(t)) for t in trajectories]
    frames = [md.load_frame(get_fn(f), index=r, top=get_fn('native.pdb')) for f, r in zip(files, rand)]

    for traj, frame, r, f in zip(trajectories, frames, rand, files):
        def test():
            eq(traj[r].xyz, frame.xyz)
            eq(traj[r].unitcell_vectors, frame.unitcell_vectors)
            eq(traj[r].time, frame.time, err_msg='%d, %d: %s' % (traj[r].time[0], frame.time[0], f))
        test.description = 'test_load_frame: %s' % f
        yield test

    t1 = md.load(get_fn('2EQQ.pdb'))
    r = np.random.randint(len(t1))
    t2 = md.load_frame(get_fn('2EQQ.pdb'), r)
    eq(t1[r].xyz, t2.xyz)


def test_iterload():
    t_ref = md.load(get_fn('frame0.h5'))[:20]
    with enter_temp_directory():
        for ext in t_ref._savers().keys():

            # only a 1 frame per file format
            if ext in ('.ncrst', '.rst7'):
                continue

            fn = 'temp%s' % ext
            t_ref.save(fn)

            def test():
                for stride in [1, 2, 3]:
                    loaded = md.load(fn, top=t_ref, stride=stride)
                    iterloaded = functools.reduce(lambda a, b: a.join(b), md.iterload(fn, top=t_ref, stride=stride, chunk=6))
                    eq(loaded.xyz, iterloaded.xyz)
                    eq(loaded.time, iterloaded.time)
                    eq(loaded.unitcell_angles, iterloaded.unitcell_angles)
                    eq(loaded.unitcell_lengths, iterloaded.unitcell_lengths)

            test.description = 'test_iterload: %s' % ext
            yield test


def test_save_load():
    # this cycles all the known formats you can save to, and then tries
    # to reload, using just a single-frame file.

    t_ref = md.load(get_fn('native.pdb'))
    t_ref.unitcell_vectors = np.array([[[1,0,0], [0,1,0], [0,0,1]]])
    with enter_temp_directory():
        for ext in t_ref._savers().keys():
            def test():
                fn = 'temp%s' % ext
                t_ref.save(fn)
                t = md.load(fn, top=t_ref.topology)

                eq(t.xyz, t_ref.xyz)
                eq(t.time, t_ref.time)
                if t._have_unitcell:
                    eq(t.unitcell_angles, t_ref.unitcell_angles)
                    eq(t.unitcell_lengths, t_ref.unitcell_lengths)

            test.description = 'test_save_load: %s' % ext
            yield test


def test_length():
    files = ['frame0.nc', 'frame0.h5', 'frame0.xtc', 'frame0.trr',
             'frame0.mdcrd', '4waters.arc', 'frame0.dcd', '2EQQ.pdb',
             'frame0.binpos', 'legacy_msmbuilder_trj0.lh5',
             'frame0.lammpstrj', 'frame0.xyz']
    for file in files:
        if file.endswith('.mdcrd'):
            kwargs = {'n_atoms': 22}
        else:
            kwargs = {}
        def f():
            try:
                eq(len(md.open(get_fn(file), **kwargs)),
                   len(md.load(get_fn(file), top=get_fn('native.pdb'))))
            except NotImplementedError as e:
                raise SkipTest(e)
        f.description = 'Length of file object: %s' % file
        yield f


def test_unitcell():
    # make sure that bogus unitcell vecotrs are not saved
    top = md.load(get_fn('native.pdb')).restrict_atoms(range(5)).topology
    t = md.Trajectory(xyz=np.random.randn(100, 5, 3), topology=top)

    for e in [tmpfns['xtc'], tmpfns['dcd'], tmpfns['binpos'], tmpfns['trr'], tmpfns['h5'], tmpfns['pdb'], tmpfns['pdb.gz'], tmpfns['nc']]:
        t.save(e)
        f = lambda: eq(md.load(e, top=top).unitcell_vectors, None)
        f.description = 'unitcell preservation in %s' % os.path.splitext(fn)[1]
        yield f

def test_chunk0_iterload():
    filename = 'frame0.h5'

    trj0 = md.load(get_fn(filename))

    for trj in md.iterload(get_fn(filename), chunk=0):
        pass

    eq(trj0.n_frames, trj.n_frames)
