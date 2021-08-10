##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christoph Klein
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

from __future__ import print_function

import numpy as np

import mdtraj as md
import mdtraj.geometry
from mdtraj.testing import eq
import pytest

"""
Reference values taken from gromacs 5.0.4 using the following commands on the
`tip3_300K_1ATM.xtc` trajectory.

echo -e "a O \n a H1 | a H2 \n q" | gmx make_ndx -f tip3p_300K_1ATM.pdb \
        -o tip3p_300K_1ATM.ndx

echo 3 3 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3_300K_1ATM_O-O_rdf.xvg
echo 3 4 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3p_300K_1ATM_O-H_rdf.xvg
echo 0 0 | gmx rdf -f tip3p_300K_1ATM.xtc -s tip3p_300K_1ATM.pdb -bin 0.005 \
        -rdf atom -n tip3p_300K_1ATM.ndx -o tip3p_300K_1ATM_all-all_rdf.xvg

"""


def test_rdf_norm(get_fn):
    # Check if the RDF's tail is normalized to ~1.0
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))

    pairs = traj.top.select_pairs('all', 'all')
    _, rdf_all = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    assert eq(np.ones(20), rdf_all[-20:], decimal=1)

    pairs = traj.top.select_pairs('name O', 'name O')
    _, rdf_O_O = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    assert eq(np.ones(20), rdf_O_O[-20:], decimal=1)

    pairs = traj.top.select_pairs('name O', "name =~ 'H.*'")
    _, rdf_O_H = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    assert eq(np.ones(20), rdf_O_H[-20:], decimal=1)


@pytest.mark.skip('Binning does not match up with gromacs currently.')
def test_compare_rdf_to_gromacs(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
    pairs = traj.top.select_pairs('all', 'all')
    R, rdf_all = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_all-all_rdf.xvg'), R, rdf_all)

    pairs = traj.top.select_pairs('name O', 'name O')
    _, rdf_O_O = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_O-O_rdf.xvg'), R, rdf_O_O)

    pairs = traj.top.select_pairs('name O', "name =~ 'H.*'")
    _, rdf_O_H = mdtraj.geometry.rdf.compute_rdf(traj, pairs)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_O-H_rdf.xvg'), R, rdf_O_H)


def compare_gromacs_xvg(ref_filename, r, g_r):
    """Compare RDFs to a reference file produced by gromacs. """
    data = np.loadtxt(ref_filename, skiprows=16)
    r0 = data[:, 0]
    g_r0 = data[:, 1]
    eq(r, r0, decimal=2)
    eq(g_r, g_r0, decimal=2)

@pytest.mark.parametrize('periodic, opt', [(True, True), (True, False), (False, True), (False, False)])
def test_rdf_t(get_fn, periodic, opt):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
    pairs = traj.top.select_pairs('name O', 'name O')
    times = list()
    for j in range(2):
        for i in range(4):
            times.append([j*4+i, j*4])

    r, g_r_t = md.compute_rdf_t(traj, pairs, times, r_range=(0, 1), periodic=periodic, opt=opt)

    assert g_r_t.shape == (8, 200)

@pytest.mark.parametrize('time_diff', [0, 2, 4])
def test_rdf_t_norm(get_fn, time_diff):
    # Check if the tail of g(r,t) at t=0 and t=time_diff is normalized to ~1.0
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))

    times = list()
    for j in range(100):
        if (j + time_diff) > traj.n_frames:
            continue
        times.append([j, j+time_diff])

    pairs = traj.top.select_pairs('name O', 'name O')
    r_t, rdf_O_O = mdtraj.geometry.rdf.compute_rdf_t(traj, pairs, times)
    assert eq(np.ones(20), np.mean(rdf_O_O, axis=0)[-20:], decimal=1)

    pairs = traj.top.select_pairs('name O', "name =~ 'H.*'")
    r_t, rdf_O_H = mdtraj.geometry.rdf.compute_rdf_t(traj, pairs, times)
    assert eq(np.ones(20), np.mean(rdf_O_H, axis=0)[-20:], decimal=1)

@pytest.mark.parametrize('periodic, opt', [(True, True), (True, False), (False, True), (False, False)])
def test_compare_rdf_t_at_0(get_fn, periodic, opt):
    # Compute g(r, t) at t = 0 and compare to g(r)
    frames = 2
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
    pairs = traj.top.select_pairs('name O', 'name O')
    times = list()
    for j in range(frames):
        times.append([j, j])

    r_t, g_r_t = md.compute_rdf_t(traj, pairs, times, self_correlation=False, periodic=periodic, opt=opt)
    mean_g_r_t = np.mean(g_r_t, axis=0)
    r, g_r = md.compute_rdf(traj[:frames], pairs, periodic=periodic, opt=periodic)

    assert eq(r_t, r)
    assert eq(mean_g_r_t, g_r, decimal=1)


@pytest.mark.parametrize('periodic, opt', [(True, True), (True, False), (False, True), (False, False)])
def test_amber_rdf_t(get_fn, periodic, opt):
    # Test triclinic case where simple approach in Tuckerman text does not
    # always work
    ext_ref = np.array([17.4835, 22.2418, 24.2910, 22.5505, 12.8686, 22.1090,
                        7.4472, 22.4253, 19.8283, 20.6935]) / 10
    traj = md.load(get_fn('test_good.nc'), top=get_fn('test.parm7'))
    pairs = [[0, 9999]]
    time_pairs = [[0, 2]]
    r_t, g_r_t = md.compute_rdf_t(traj, pairs, time_pairs, 
            self_correlation=False, periodic=periodic, opt=opt)


@pytest.mark.skip('Binning does not match up with gromacs currently.')
def test_compare_gromacs_rdf_t(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))
    pairs = traj.top.select_pairs('name O', 'name O')
    times = list()
    for j in range(traj.n_frames):
        times.append([j, j])

    r_t, g_r_t = md.compute_rdf_t(traj, pairs, times, self_correlation=False, periodic=True, opt=True)
    mean_g_r_t = np.mean(g_r_t, axis=0)
    compare_gromacs_xvg(get_fn('tip3p_300K_1ATM_O-O_rdf.xvg'), r_t, mean_g_r_t)


def test_compare_rdf_t_master(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))

    times = list()
    for j in range(100):
        times.append([0, j])

    pairs = traj.top.select_pairs('name O', 'name O')
    r_t, rdf_O_O = mdtraj.geometry.rdf.compute_rdf_t(traj, pairs, times)

    master_r_t = np.loadtxt(get_fn('r_O_O_rdf_t.txt'))
    master_g_r_t = np.loadtxt(get_fn('O_O_rdf_t.txt'))

    assert eq(r_t, master_r_t)
    assert eq(rdf_O_O, master_g_r_t, decimal=5)

def test_compare_n_concurrent_pairs(get_fn):
    traj = md.load(get_fn('tip3p_300K_1ATM.xtc'), top=get_fn('tip3p_300K_1ATM.pdb'))

    times = list()
    for j in range(100):
        times.append([0, j])

    pairs = traj.top.select_pairs('name O', 'name O')

    r_t_1, rdf_O_O_1 = mdtraj.geometry.rdf.compute_rdf_t(traj, pairs, times, n_concurrent_pairs=100000)
    r_t_2, rdf_O_O_2 = mdtraj.geometry.rdf.compute_rdf_t(traj, pairs, times, n_concurrent_pairs=314159)

    assert eq(r_t_1, r_t_2)
    assert eq(rdf_O_O_1, rdf_O_O_2)
