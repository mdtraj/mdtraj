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


##############################################################################
# Imports
##############################################################################

import numpy as np
from numpy.testing import *

import mdtraj as md
from mdtraj.testing import get_fn
from mdtraj.geometry.sasa import _ATOMIC_RADII

##############################################################################
# Globals
##############################################################################


# set up a mock topology with 1 atom
topology1 = md.Topology()
topology1.add_atom('H', md.pdb.element.hydrogen, topology1.add_residue('res', topology1.add_chain()))

# set up a mock topology with two atoms
topology2 = md.Topology()
_res2 = topology2.add_residue('res', topology2.add_chain())
topology2.add_atom('H', md.pdb.element.hydrogen, _res2)
topology2.add_atom('H', md.pdb.element.hydrogen, _res2)


##############################################################################
# Tests
##############################################################################


def test_sasa_0():
    # make one atom at the origin
    traj = md.Trajectory(xyz=np.zeros((1,1,3)), topology=topology1)

    probe_radius = 0.14
    calc_area = np.sum(md.geometry.shrake_rupley(traj, probe_radius=probe_radius))
    true_area = 4 * np.pi * (_ATOMIC_RADII['H'] + probe_radius)**2

    assert_approx_equal(calc_area, true_area)


def test_sasa_1():
    # two atoms
    traj = md.Trajectory(xyz=np.zeros((1,2,3)), topology=topology2)

    probe_radius = 0.14
    true  = 4 * np.pi * (_ATOMIC_RADII['H'] + probe_radius)**2

    # when atoms are closer than 2e-5, there seems to be a bug.
    # note that you should never actually have a case where atoms are this close
    # but nonetheless I'm adding a check for this in the implementation -- to make
    # it crash if the atoms are too close, as opposed to giving you wrong results
    separations = np.linspace(2.0e-5, probe_radius*2 + _ATOMIC_RADII['H']*2, 10)
    areas = np.zeros_like(separations)

    # check the sasa as we vary the separation
    for i, sep in enumerate(separations):
        traj.xyz[0, 0, 1] = sep
        areas[i] = np.sum(md.geometry.shrake_rupley(traj, probe_radius=probe_radius))

    assert_approx_equal(areas[0], true)
    assert_approx_equal(areas[-1], 2*true)
    # make sure that areas is increasing
    assert_array_less(areas[0:8], areas[1:9])


def test_sasa_2():
    t = md.load(get_fn('frame0.h5'))
    val1 = np.sum(md.geometry.shrake_rupley(t[0])) # calculate only frame 0
    val2 = np.sum(md.geometry.shrake_rupley(t)[0]) # calculate on all frames
    true_frame_0_sasa = 2.859646797180176

    assert_approx_equal(true_frame_0_sasa, val1)
    assert_approx_equal(true_frame_0_sasa, val2)


def test_sasa_3():
    traj_ref = np.loadtxt(get_fn('g_sas_ref.dat'))
    traj = md.load(get_fn('frame0.h5'))
    traj_sasa = md.geometry.shrake_rupley(traj, probe_radius=0.14, n_sphere_points = 960)

    # the algorithm used by gromacs' g_sas is slightly different than the one
    # used here, so the results are not exactly the same
    assert_array_almost_equal(traj_sasa, traj_ref, decimal=2)
