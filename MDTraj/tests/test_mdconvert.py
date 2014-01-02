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


"""
Code to test the mdconvert script. These tests take about two minutes to run.

This checks all pairs for formats, converting from format x -> format y. it
also trys using striding to subsample the trajectory and atom_indices, so it
does significant integration testing of the XXXTrajectoryFile modules as well.
"""

##############################################################################
# imports
##############################################################################

import os
import tempfile
import shutil
import numpy as np

import mdtraj as md
from mdtraj.utils import import_
from mdtraj.testing import skipif, get_fn, eq, slow

try:
    scripttest = import_('scripttest')
    HAVE_SCRIPTTEST = True
except SystemExit:
    HAVE_SCRIPTTEST = False

##############################################################################
# globals
##############################################################################

# if you switch DEBUG_MODE to True, none of the files will deleted
# at the end of the execution of this suite, so that you can debug the
# problem by running mdconvert manually.
DEBUG_MODE = False
# DEBUG_MODE = False

staging_dir = tempfile.mkdtemp()
output_dir = os.path.join(staging_dir, 'output')
def teardown_module(module):
    if not DEBUG_MODE:
        shutil.rmtree(staging_dir)

def setup_module():
    global TRAJ

    xyz = np.around(np.random.randn(10, 5, 3).astype(np.float32), 2)
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue('ALA', chain)
    topology.add_atom('CA', md.pdb.element.carbon, residue)
    topology.add_atom('HG1', md.pdb.element.hydrogen, residue)
    topology.add_atom('SG', md.pdb.element.sulfur, residue)
    topology.add_atom('OD1', md.pdb.element.oxygen, residue)
    topology.add_atom('NE', md.pdb.element.nitrogen, residue)

    time = np.arange(10)**2
    unitcell_lengths = np.array([[1.1,1.2,1.3]] * 10)
    unitcell_angles = np.array([[90, 90, 95]] * 10)

    TRAJ = md.Trajectory(xyz, topology=topology, time=time,
                         unitcell_lengths=unitcell_lengths,
                         unitcell_angles=unitcell_angles)


##############################################################################
# test
##############################################################################


@skipif(not HAVE_SCRIPTTEST)
def test_mdconvert_index():
    "Check that extracting a specific index works"
    env = scripttest.TestFileEnvironment(output_dir)
    path = os.path.join(staging_dir, 'traj.h5')
    TRAJ.save(path)
    command = ['mdconvert', path, '-i 4', '-o', 'frame4.pdb']
    env.run(*command, expect_stderr=True)
    frame4 = md.load(os.path.join(output_dir, 'frame4.pdb'))
    eq(frame4.xyz, TRAJ[4].xyz)

    os.unlink(path)


@skipif(not HAVE_SCRIPTTEST)
def test_mdconvert_slice():
    "Check that extracting a specific slice works"
    env = scripttest.TestFileEnvironment(output_dir)
    path = os.path.join(staging_dir, 'traj.h5')
    TRAJ.save(path)
    command = ['mdconvert', path, '-i 1:5:2', '-o', 'frame13.pdb']
    env.run(*command, expect_stderr=True)
    frame13 = md.load(os.path.join(output_dir, 'frame13.pdb'))
    eq(frame13.xyz, TRAJ[1:5:2].xyz)
    os.unlink(path)


@slow
@skipif(not HAVE_SCRIPTTEST)
def test_mdconvert_0():
    """ensure that the xyz coordinates are preserved by a trip
       from python -> save in format X -> mdconvert to format Y -> python
    """
    env = scripttest.TestFileEnvironment(output_dir)

    # save one copy of traj for use as a topology file
    topology_fn = os.path.join(staging_dir, 'topology.pdb')
    TRAJ[0].save(topology_fn)

    # save a .dat file for the atom_indices so that we can test
    # mdconvert's atom_indices flag
    atom_indices = np.array([0, 3])
    atom_indices_fn = os.path.join(staging_dir, 'atom_indices.dat')
    np.savetxt(atom_indices_fn, atom_indices, fmt='%d')

    # fns = ['traj.xtc', 'traj.dcd', 'traj.binpos', 'traj.trr', 'traj.nc',
    #        'traj.pdb', 'traj.h5', 'traj.lh5']
    fns = ['traj.xtc', 'traj.lh5']

    for fn in fns:
        path = os.path.join(staging_dir, fn)
        TRAJ.save(path)

        for fn2 in filter(lambda e: e != fn, fns):
            ext1, ext2 = [os.path.splitext(f)[1] for f in [fn, fn2]]

            command1 = ['mdconvert', path, '-o', fn2, '-c 6']
            if ext2 in ['.pdb', '.h5', '.lh5']:
                # if we're saving a pdb or h5, we need to give it a topology too
                command1 += ['-t', topology_fn]

            # one set of tests, with no extra flags to mdconvert
            execution1 = lambda : env.run(*command1, expect_stderr=True)
            execution1.description = 'mdconvert: converting %s -> %s' % (fn, fn2)

            # lets try using the --atom_indices flag to mdconvert
            command2 = command1 + ['-a', atom_indices_fn]
            command2[3] = 'subset.' + fn2   # make sure the output goes to a different file
            execution2 = lambda : env.run(*command2, expect_stderr=True)
            execution2.description = 'mdconvert: converting %s -> %s (atom_indices)' % (fn, 'subset.' + fn2)

            # lets try one using the --stride 3 flag
            command3 = command1 +  ['-s 3']
            command3[3] = 'stride.' + fn2   # change the out filename, so they don't clobbed
            execution3 = lambda : env.run(*command3, expect_stderr=True)
            execution3.description = 'mdconvert: converting %s -> %s (stride)' % (fn, 'stride.' + fn2)

            yield execution1
            yield execution2
            yield execution3

            # ensure that the xyz coordinates are preserved by a trip
            # from python -> save in format X -> mdconvert to format Y -> python
            load_kwargs_check1, load_kwargs_check2 = {}, {}
            if ext2 not in ['.pdb', '.h5', '.lh5']:
                load_kwargs_check1['top'] = TRAJ.topology
                load_kwargs_check2['top'] = TRAJ.topology.subset(atom_indices)

            def check():
                out1 = md.load(os.path.join(output_dir, fn2), **load_kwargs_check1)
                out2 = md.load(os.path.join(output_dir, 'subset.' + fn2), **load_kwargs_check2)
                out3 = md.load(os.path.join(output_dir, 'stride.' + fn2), **load_kwargs_check1)

                eq(out1.xyz, TRAJ.xyz)
                eq(out2.xyz, TRAJ.xyz[:, atom_indices])
                eq(out3.xyz, TRAJ.xyz[::3])

                if ext1 not in ['.binpos', '.lh5'] and ext2 not in ['.binpos', '.lh5']:
                    # binpos doesn't save unitcell information
                    eq(out1.unitcell_vectors, TRAJ.unitcell_vectors, decimal=2)
                    eq(out2.unitcell_vectors, TRAJ.unitcell_vectors, decimal=2)
                    eq(out3.unitcell_vectors, TRAJ.unitcell_vectors[::3], decimal=2)

                if all(e in ['.xtc', '.trr', '.nc', '.h5'] for e in [ext1, ext2]):
                    # these formats contain time information
                    eq(out1.time, TRAJ.time)
                    eq(out2.time, TRAJ.time)
                    eq(out3.time, TRAJ.time[::3])

                if ext2 in ['.pdb', '.h5', '.lh5']:
                    # these formats contain a topology in the file that was
                    # read from disk
                    eq(out1.topology, TRAJ.topology)
                    eq(out2.topology, TRAJ.topology.subset(atom_indices))
                    eq(out3.topology, TRAJ.topology)

            check.description = 'mdconvert: checking %s -> %s' % (fn, fn2)
            yield check

            if not DEBUG_MODE:
                os.unlink(os.path.join(output_dir, fn2))
                os.unlink(os.path.join(output_dir, 'subset.' + fn2))
                os.unlink(os.path.join(output_dir, 'stride.' + fn2))

        if not DEBUG_MODE:
            os.unlink(path)
