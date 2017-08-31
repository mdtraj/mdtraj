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
also tries using striding to subsample the trajectory and atom_indices, so it
does significant integration testing of the XXXTrajectoryFile modules as well.
"""

import os
import sys
import numpy as np

import mdtraj as md
from mdtraj import element
from mdtraj.testing import eq
import pytest
import subprocess

on_win = (sys.platform == 'win32')
on_py3 = (sys.version_info >= (3, 0))


@pytest.fixture()
def traj(tmpdir):
    xyz = np.around(np.random.randn(10, 5, 3).astype(np.float32), 2)
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue('ALA', chain)
    topology.add_atom('CA', element.carbon, residue)
    topology.add_atom('HG1', element.hydrogen, residue)
    topology.add_atom('SG', element.sulfur, residue)
    topology.add_atom('OD1', element.oxygen, residue)
    topology.add_atom('NE', element.nitrogen, residue)

    time = np.arange(10) ** 2
    unitcell_lengths = np.array([[1.1, 1.2, 1.3]] * 10)
    unitcell_angles = np.array([[90, 90, 95]] * 10)

    traj = md.Trajectory(xyz, topology=topology, time=time,
                         unitcell_lengths=unitcell_lengths,
                         unitcell_angles=unitcell_angles)

    fn = '{}/ref.h5'.format(tmpdir)
    traj.save(fn)
    return traj, fn, str(tmpdir)


def test_index(traj):
    # Check that extracting a specific index works
    traj, in_fn, tmpdir = traj
    out_fn = '{}/frame4.pdb'.format(tmpdir)
    subprocess.check_call(['mdconvert', in_fn, '-i', '4', '-o', out_fn])
    frame4 = md.load(out_fn)
    eq(frame4.xyz, traj[4].xyz)


def test_slice(traj):
    # Check that extracting a specific slice works
    traj, in_fn, tmpdir = traj
    out_fn = '{}/frame13.pdb'.format(tmpdir)
    subprocess.check_call(['mdconvert', in_fn, '-i', '1:5:2', '-o', out_fn])
    frame13 = md.load(out_fn)
    eq(frame13.xyz, traj[1:5:2].xyz)


extensions = [
    'xtc', 'dcd', 'binpos', 'trr', 'nc',
    'pdb', 'h5', 'lh5', 'netcdf'
]


@pytest.fixture(params=extensions, ids=lambda x: 'from-' + x)
def extension(request):
    if on_win and on_py3 and request.param == 'lh5':
        pytest.skip('No lh5 on windows py3')
    return request.param


def test_pairwise(traj, extension):
    """ensure that the xyz coordinates are preserved by a trip
       from python -> save in format X -> mdconvert to format Y -> python
    """
    traj, _, tmpdir = traj
    ext1 = extension

    # save one copy of traj for use as a topology file
    topology_fn = "{}/topology.pdb".format(tmpdir)
    traj[0].save(topology_fn)

    # save a .dat file for the atom_indices so that we can test
    # mdconvert's atom_indices flag
    atom_indices = np.array([0, 3])
    atom_indices_fn = "{}/atom_indices.dat".format(tmpdir)
    np.savetxt(atom_indices_fn, atom_indices, fmt='%d')

    in_fn = "{}/traj.{}".format(tmpdir, ext1)
    traj.save(in_fn)
    working_dir = '{}/from-{}'.format(tmpdir, ext1)
    os.mkdir(working_dir)

    for ext2 in extensions:
        print(ext2)
        out_fn = 'traj.{}'.format(ext2)

        command1 = ['mdconvert', in_fn, '-o', out_fn, '-c 6']
        if ext2 in ['pdb', 'h5', 'lh5']:
            # if we're saving a pdb or h5, we need to give it a topology too
            command1 += ['-t', topology_fn]

        # TODO: test fixture
        subprocess.check_call(command1, cwd=working_dir)

        # Use the --atom_indices flag to mdconvert
        command2 = command1 + ['-a', atom_indices_fn]
        command2[3] = 'subset.' + out_fn  # make sure the output goes to a different file
        subprocess.check_call(command2, cwd=working_dir)

        # Use the --stride 3 flag
        command3 = command1 + ['-s 3']
        command3[3] = 'stride.' + out_fn  # change the out filename, so they don't clobbed
        subprocess.check_call(command3, cwd=working_dir)

        # ensure that the xyz coordinates are preserved by a trip
        # from python -> save in format X -> mdconvert to format Y -> python
        load_kwargs_check1 = {}
        load_kwargs_check2 = {}
        if ext2 not in ['pdb', 'h5', 'lh5']:
            load_kwargs_check1['top'] = traj.topology
            load_kwargs_check2['top'] = traj.topology.subset(atom_indices)

        out1 = md.load(os.path.join(working_dir, out_fn), **load_kwargs_check1)
        out2 = md.load(os.path.join(working_dir, 'subset.' + out_fn), **load_kwargs_check2)
        out3 = md.load(os.path.join(working_dir, 'stride.' + out_fn), **load_kwargs_check1)

        if ext1 in ['lh5'] or ext2 in ['lh5']:
            decimal = 3
        else:
            decimal = 6
        eq(out1.xyz, traj.xyz, decimal=decimal)
        eq(out2.xyz, traj.xyz[:, atom_indices], decimal=decimal)
        eq(out3.xyz, traj.xyz[::3], decimal=decimal)

        if ext1 not in ['binpos', 'lh5'] and ext2 not in ['binpos', 'lh5']:
            # binpos doesn't save unitcell information
            eq(out1.unitcell_vectors, traj.unitcell_vectors, decimal=2)
            eq(out2.unitcell_vectors, traj.unitcell_vectors, decimal=2)
            eq(out3.unitcell_vectors, traj.unitcell_vectors[::3], decimal=2)

        if all(e in ['xtc', 'trr', 'nc', 'h5'] for e in [ext1, ext2]):
            # these formats contain time information
            eq(out1.time, traj.time)
            eq(out2.time, traj.time)
            eq(out3.time, traj.time[::3])

        if ext2 in ['pdb', 'h5', 'lh5']:
            # these formats contain a topology in the file that was
            # read from disk
            eq(out1.topology, traj.topology)
            eq(out2.topology, traj.topology.subset(atom_indices))
            eq(out3.topology, traj.topology)


def test_mdconvert_alanine(tmpdir, get_fn):
    command = ['mdconvert', get_fn('alanine-dipeptide-explicit.binpos'),
               '--top', get_fn('alanine-dipeptide-explicit.prmtop'),
               '-o', 'out.dcd']
    subprocess.check_call(command, cwd=str(tmpdir))
    t = md.load('{}/out.dcd'.format(tmpdir), top=get_fn('alanine-dipeptide-explicit.prmtop'))
    t2 = md.load(get_fn('alanine-dipeptide-explicit.binpos'), top=get_fn('alanine-dipeptide-explicit.prmtop'))

    eq(t.xyz, t2.xyz)
    eq(t.topology, t2.topology)
