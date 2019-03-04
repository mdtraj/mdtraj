##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
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


import mdtraj as md
import numpy as np
from mdtraj import element
from mdtraj.geometry.sasa import _ATOMIC_RADII
from mdtraj.testing import eq

# set up a mock topology with 1 atom
topology1 = md.Topology()
topology1.add_atom('H', element.hydrogen, topology1.add_residue('res', topology1.add_chain()))

# set up a mock topology with two atoms
topology2 = md.Topology()
_res2 = topology2.add_residue('res', topology2.add_chain())
topology2.add_atom('H', element.hydrogen, _res2)
topology2.add_atom('H', element.hydrogen, _res2)

# set up a mock topologies with a single chlorine and sulfur atom
topology3 = md.Topology()
topology4 = md.Topology()
topology3.add_atom('Cl', element.chlorine, topology3.add_residue('res', topology3.add_chain()))
topology4.add_atom('S', element.sulfur, topology4.add_residue('res', topology4.add_chain()))


def test_sasa_0():
    # make one atom at the origin
    traj = md.Trajectory(xyz=np.zeros((1, 1, 3)), topology=topology1)

    probe_radius = 0.14
    calc_area = np.sum(md.geometry.shrake_rupley(traj, probe_radius=probe_radius))
    true_area = 4 * np.pi * (_ATOMIC_RADII['H'] + probe_radius) ** 2

    np.testing.assert_approx_equal(calc_area, true_area)


def test_sasa_1():
    # two atoms
    traj = md.Trajectory(xyz=np.zeros((1, 2, 3)), topology=topology2)

    probe_radius = 0.14
    true = 4 * np.pi * (_ATOMIC_RADII['H'] + probe_radius) ** 2

    # when atoms are closer than 2e-5, there seems to be a bug.
    # note that you should never actually have a case where atoms are this close
    # but nonetheless I'm adding a check for this in the implementation -- to make
    # it crash if the atoms are too close, as opposed to giving you wrong results
    separations = np.linspace(2.0e-5, probe_radius * 2 + _ATOMIC_RADII['H'] * 2, 10)
    areas = np.zeros_like(separations)

    # check the sasa as we vary the separation
    for i, sep in enumerate(separations):
        traj.xyz[0, 0, 1] = sep
        areas[i] = np.sum(md.geometry.shrake_rupley(traj, probe_radius=probe_radius))

    np.testing.assert_approx_equal(areas[0], true, significant=3)
    np.testing.assert_approx_equal(areas[-1], 2 * true)
    # make sure that areas is increasing
    np.testing.assert_array_less(areas[0:8], areas[1:9])


def test_sasa_2(get_fn):
    t = md.load(get_fn('frame0.h5'))
    val1 = np.sum(md.geometry.shrake_rupley(t[0]))  # calculate only frame 0
    val2 = np.sum(md.geometry.shrake_rupley(t)[0])  # calculate on all frames
    true_frame_0_sasa = 3.5556480884552

    np.testing.assert_approx_equal(true_frame_0_sasa, val1)
    np.testing.assert_approx_equal(true_frame_0_sasa, val2)


def test_sasa_3(get_fn):
    traj_ref = np.loadtxt(get_fn('gmx_sasa.dat'))
    traj = md.load(get_fn('frame0.h5'))
    traj_sasa = md.geometry.shrake_rupley(traj, probe_radius=0.14, n_sphere_points=960)

    # the algorithm used by gromacs' g_sas is slightly different than the one
    # used here, so the results are not exactly the same
    np.testing.assert_array_almost_equal(traj_sasa, traj_ref, decimal=1)


def test_sasa_4(get_fn):
    def _test_atom_group(t, value):
        sasa = md.shrake_rupley(t, mode='atom')
        rids = np.array([a.residue.index for a in t.top.atoms])

        for i, rid in enumerate(np.unique(rids)):
            mask = (rids == rid)
            eq(value[:, i], np.sum(sasa[:, mask], axis=1))

    t = md.load(get_fn('frame0.h5'))
    value = md.shrake_rupley(t, mode='residue')
    assert value.shape == (t.n_frames, t.n_residues)
    _test_atom_group(t, value)

    # scramle the order of the atoms, and which residue each is a
    # member of
    df, bonds = t.top.to_dataframe()
    df['resSeq'] = np.random.permutation(df['resSeq'])
    df['resName'] = df['resSeq']
    t.top = md.Topology.from_dataframe(df, bonds)

    value = md.shrake_rupley(t, mode='residue')
    _test_atom_group(t, value)


def test_sasa_5():
    # Test if you can change atom radii effectively without breaking things
    traj_chlorine = md.Trajectory(xyz=np.zeros((1, 1, 3)), topology=topology3)
    traj_sulfur = md.Trajectory(xyz=np.zeros((1, 1, 3)), topology=topology4)

    probe_radius = 0.14
    change_radii = {'Cl': 0.175}
    calc_area_chlorine = np.sum(md.geometry.shrake_rupley(traj_chlorine,
                                                          probe_radius=probe_radius,
                                                          change_radii=change_radii))
    calc_area_sulfur = np.sum(md.geometry.shrake_rupley(traj_sulfur,
                                                        probe_radius=probe_radius,
                                                        change_radii=change_radii))
    true_area_chlorine = 4 * np.pi * (0.175 + probe_radius) ** 2
    # This atom should remain unchanged
    true_area_sulfur = 4 * np.pi * (_ATOMIC_RADII['S'] + probe_radius) ** 2
    # And calc_area_chlorine should be different from this value
    prev_area_chlorine = 4 * np.pi * (_ATOMIC_RADII['Cl'] + probe_radius) ** 2

    np.testing.assert_approx_equal(calc_area_chlorine, true_area_chlorine)
    # Also ensure that rest of dict is not gettting deleted
    np.testing.assert_approx_equal(calc_area_sulfur, true_area_sulfur)
    # Finally, ensure that we are not changing _ATOMIC_RADII
    assert prev_area_chlorine != calc_area_chlorine
