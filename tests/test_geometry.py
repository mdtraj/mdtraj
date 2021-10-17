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


import itertools

import mdtraj as md
import mdtraj.geometry
import numpy as np
import pytest
from mdtraj.testing import eq

RUN_PERFORMANCE_TESTS = False


def test_rg(get_fn):
    t0 = md.load(get_fn('traj.h5'))
    Rg = mdtraj.geometry.rg.compute_rg(t0)
    Rg0 = np.loadtxt(get_fn("Rg_traj_ref.dat"))
    eq(Rg, Rg0)


def test_distances(get_fn):
    t0 = md.load(get_fn('traj.h5'))
    atom_pairs = np.loadtxt(get_fn("atom_pairs.dat"), 'int')
    distances = mdtraj.geometry.distance.compute_distances(t0, atom_pairs)
    distances0 = np.loadtxt(get_fn("atom_distances_traj_ref.dat"))
    eq(distances, distances0)


def test_center_of_mass():
    top = md.Topology()
    chain = top.add_chain()
    resi = top.add_residue("NA", chain)
    top.add_atom('H1', md.element.hydrogen, resi)
    top.add_atom('H2', md.element.hydrogen, resi)
    top.add_atom('O', md.element.oxygen, resi)

    xyz = np.array([
        [
            # Frame 1 - Right triangle
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]
        ],
        [
            # Frame 2
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.5, 0.0]
        ]
    ])
    traj = md.Trajectory(xyz, top)
    com_test = mdtraj.geometry.distance.compute_center_of_mass(traj)

    com_actual = (1 / 18.015324) * np.array([
        [1.007947, 1.007947, 0.0],
        [1.007947 + 0.5 * 15.99943, 1.007947 + 0.5 * 15.99943, 0.0],
    ])

    eq(com_actual, com_test, decimal=4)

    #test with selection string
    com_test = mdtraj.geometry.distance.compute_center_of_mass(traj, select='index 0')

    com_actual = np.array([
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ])

    eq(com_actual, com_test, decimal=4)


def test_center_of_geometry():
    top = md.Topology()
    chain = top.add_chain()
    resi = top.add_residue("NA", chain)
    top.add_atom('H1', md.element.hydrogen, resi)
    top.add_atom('H2', md.element.hydrogen, resi)
    top.add_atom('O', md.element.oxygen, resi)

    xyz = np.array([
        [
            # Frame 1 - Right triangle
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]
        ],
        [
            # Frame 2
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.5, 0.0]
        ]
    ])
    traj = md.Trajectory(xyz, top)
    center_test = mdtraj.geometry.distance.compute_center_of_geometry(traj)

    center_actual = np.array([
        [1.0 / 3, 1.0 / 3, 0.0],
        [1.5 / 3, 1.5 / 3, 0.0],
    ])

    eq(center_actual, center_test, decimal=4)


def test_dihedral_indices(get_fn):
    traj = md.load(get_fn('1bpi.pdb'))
    # Manually compare generated indices to known answers.
    phi0_ind = np.array([3, 12, 13, 14]) - 1  # Taken from PDB, so subtract 1
    psi0_ind = np.array([1, 2, 3, 12]) - 1  # Taken from PDB, so subtract 1

    ind = mdtraj.geometry.dihedral.indices_phi(traj)
    eq(ind[0], phi0_ind)

    ind = mdtraj.geometry.dihedral.indices_psi(traj)
    eq(ind[0], psi0_ind)


def test_dihedral_index_offset_generation(get_fn):
    traj = md.load(get_fn('1bpi.pdb'))
    top = traj.topology

    # The atom indices of the first phi angle
    result = np.array([2, 11, 12, 13])
    print([e.name for e in traj.topology.atoms])
    ind1 = mdtraj.geometry.dihedral.indices_phi(top)
    rid2, ind2 = mdtraj.geometry.dihedral._atom_sequence(top, ["-C", "N", "CA", "C"])
    rid3, ind3 = mdtraj.geometry.dihedral._atom_sequence(top, ["-C", "N", "CA", "C"], [-1, 0, 0, 0])
    rid4, ind4 = mdtraj.geometry.dihedral._atom_sequence(top, ["C", "N", "CA", "C"], [-1, 0, 0, 0])

    eq(ind1, ind2)
    eq(ind1, ind3)
    eq(ind1, ind4)
    eq(ind1[0], result)


def test_dihedral_0(get_fn):
    """We compared phi and psi angles from pymol to MDTraj output."""
    traj = md.load(get_fn('1bpi.pdb'))[0]
    indices, phi = mdtraj.geometry.dihedral.compute_phi(traj)
    phi0 = np.array([-34.50956, -50.869690]) * np.pi / 180.  # Pymol
    eq(phi[0, 0:2], phi0, decimal=4)
    eq(indices[0], np.array([2, 11, 12, 13]))

    indices, psi = mdtraj.geometry.dihedral.compute_psi(traj)
    psi0 = np.array([134.52554, 144.880173]) * np.pi / 180.  # Pymol
    eq(psi[0, 0:2], psi0, decimal=4)
    eq(indices[0], np.array([0, 1, 2, 11]))

    indices, chi = mdtraj.geometry.dihedral.compute_chi1(traj)
    chi0 = np.array([-43.37841, -18.14592]) * np.pi / 180.  # Pymol
    eq(chi[0, 0:2], chi0, decimal=4)
    eq(indices[0], np.array([0, 1, 4, 5]))


def test_dihedral_1():
    n_atoms = 4
    np.random.seed(42)
    xyz = np.random.randn(5, n_atoms, 3)
    t = md.Trajectory(xyz=xyz, topology=None)
    indices = list(itertools.combinations(range(n_atoms), 4))
    r1 = md.geometry.compute_dihedrals(t, indices, opt=False)
    r2 = md.geometry.compute_dihedrals(t, indices, opt=True)
    eq(r1, r2)


def test_dihedral_triclinic():
    n_atoms = 4
    np.random.seed(42)
    xyz = np.random.randn(5, n_atoms, 3)
    t = md.Trajectory(xyz=xyz, topology=None,
                      unitcell_lengths=np.outer(np.ones(5), np.array([1.5, 1.5, 1.5])),
                      unitcell_angles=np.outer(np.ones(5), np.array([95, 100, 110])))
    indices = list(itertools.combinations(range(n_atoms), 4))
    r1 = md.geometry.compute_dihedrals(t, indices, opt=False)
    r2 = md.geometry.compute_dihedrals(t, indices, opt=True)
    eq(r1, r2)


def test_angle_0():
    xyz = np.array([[[0, 0, 0],
                     [0, 1, 0],
                     [1, 1, 0]]])
    t = md.Trajectory(xyz=xyz, topology=None)
    result = np.array(np.pi / 2).reshape(1, 1)
    eq(result, md.geometry.compute_angles(t, [[0, 1, 2]], opt=False))
    eq(result, md.geometry.compute_angles(t, [[0, 1, 2]], opt=True))


def test_angle_periodic_0():
    xyz = np.array([[[0, 0, 0],
                     [0, 1, 0],
                     [1, 1, 0]]])
    t = md.Trajectory(xyz=xyz, topology=None,
                      unitcell_lengths=np.array([10, 10, 10]),
                      unitcell_angles=np.array([90, 90, 90]))
    result = np.array(np.pi / 2).reshape(1, 1)
    eq(result, md.geometry.compute_angles(t, [[0, 1, 2]], opt=False))
    eq(result, md.geometry.compute_angles(t, [[0, 1, 2]], opt=True))


def test_angle_triclinic_0():
    xyz = np.array([[[0, 0, 0],
                     [0, 1, 0],
                     [1, 1, 0]]])
    t = md.Trajectory(xyz=xyz, topology=None,
                      unitcell_lengths=np.array([1.5, 1.5, 1.5]),
                      unitcell_angles=np.array([95, 100, 110]))
    eq(md.geometry.compute_angles(t, [[0, 1, 2]], opt=False), md.geometry.compute_angles(t, [[0, 1, 2]], opt=True))


def test_angle_1():
    # the two routines sometimes give slightly different answers for
    # wierd angles really close to 0 or 180. I suspect this is because
    # different implementations of acos -- which are basically implemented
    # taylor series expansions -- are parameterized slightly differently on
    # different libraries/platforms. setting the random number seed helps to
    # ensure that this doesn't break stochastically too often.

    n_atoms = 5
    np.random.seed(24)
    xyz = np.random.randn(50, n_atoms, 3)
    t = md.Trajectory(xyz=xyz, topology=None)
    indices = list(itertools.combinations(range(n_atoms), 3))
    r1 = md.geometry.compute_angles(t, indices, opt=False)
    r2 = md.geometry.compute_angles(t, indices, opt=True)
    assert np.max(np.abs(r1 - r2)) < 1e-2
    assert np.nansum(np.abs(r1 - r2)) / r1.size < 5e-4


@pytest.mark.skipif(not RUN_PERFORMANCE_TESTS, reason='Not doing performance testing')
def test_dihedral_performance():
    n_atoms = 10
    xyz = np.random.randn(5000, n_atoms, 3)
    t = md.Trajectory(xyz=xyz, topology=None)
    indices = np.asarray(list(itertools.combinations(range(n_atoms), 4)), dtype=np.int32)
    import time
    t1 = time.time()
    r1 = md.geometry.compute_dihedrals(t, indices, opt=False)
    t2 = time.time()
    r2 = md.geometry.compute_dihedrals(t, indices, opt=True)
    t3 = time.time()

    print('\ndihedral performance:')
    print('numpy:   %f s' % (t2 - t1))
    print('opt sse: %f s' % (t3 - t2))


@pytest.mark.skipif(not RUN_PERFORMANCE_TESTS, reason='Not doing performance testing')
def test_angle_performance():
    n_atoms = 20
    xyz = np.random.randn(10000, n_atoms, 3)
    t = md.Trajectory(xyz=xyz, topology=None)
    indices = np.asarray(list(itertools.combinations(range(n_atoms), 3)), dtype=np.int32)
    import time
    t1 = time.time()
    r1 = md.geometry.compute_angles(t, indices, opt=False)
    t2 = time.time()
    r2 = md.geometry.compute_angles(t, indices, opt=True)
    t3 = time.time()

    print('\nangle performance:')
    print('numpy:   %f s' % (t2 - t1))
    print('opt sse: %f s' % (t3 - t2))


def test_angle_nan():
    t = md.Trajectory(topology=None, xyz=np.array([
        [0, 0, 0.2],
        [0, 0, 0.3],
        [0, 0, 0.0]
    ]))
    angles = md.compute_angles(t, [[0, 1, 2]], opt=True)
    np.testing.assert_array_almost_equal(angles, [[0]])


def test_angle_pbc(get_fn):
    traj_uncorrected = md.load(get_fn('1am7_uncorrected.xtc'), top=get_fn('1am7_protein.pdb'))
    traj_corrected = md.load(get_fn('1am7_corrected.xtc'), top=get_fn('1am7_protein.pdb'))

    indices = []
    for i in range(traj_uncorrected.n_residues):
        r = traj_uncorrected.topology.residue(i)
        indices.append((r.atom(0).index, r.atom(1).index, r.atom(2).index))

    epsilon = 6e-3

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=False, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=False, periodic=True)
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=True, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=True, periodic=True)
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=True, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=True, periodic=False)
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=False, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=False, periodic=False)
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=True, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=False, periodic=False)
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=False, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=True, periodic=False)
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=False, periodic=True)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=True, periodic=True)
    assert np.max(np.abs(ang1 - ang2)) < epsilon


def test_dihedral_pbc(get_fn):
    traj_uncorrected = md.load(get_fn('1am7_uncorrected.xtc'), top=get_fn('1am7_protein.pdb'))
    traj_corrected = md.load(get_fn('1am7_corrected.xtc'), top=get_fn('1am7_protein.pdb'))

    epsilon = 6E-3

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=False, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=False, periodic=True)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=True, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=True, periodic=True)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=True, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=True, periodic=False)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=False, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=False, periodic=False)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=True, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=False, periodic=False)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=False, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=True, periodic=True)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=False, periodic=True)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=True, periodic=False)[1]
    assert np.max(np.abs(ang1 - ang2)) < epsilon


def test_dihedral_pbc_fails1(get_fn):
    traj_uncorrected = md.load(get_fn('1am7_uncorrected.xtc'), top=get_fn('1am7_protein.pdb'))
    traj_corrected = md.load(get_fn('1am7_corrected.xtc'), top=get_fn('1am7_protein.pdb'))

    epsilon = 1E-2

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=False, periodic=False)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=False, periodic=False)[1]
    assert not (np.max(np.abs(ang1 - ang2)) < epsilon)


def test_dihedral_pbc_fails2(get_fn):
    traj_uncorrected = md.load(get_fn('1am7_uncorrected.xtc'), top=get_fn('1am7_protein.pdb'))
    traj_corrected = md.load(get_fn('1am7_corrected.xtc'), top=get_fn('1am7_protein.pdb'))

    epsilon = 1E-2

    ang1 = md.geometry.compute_phi(traj_uncorrected, opt=True, periodic=False)[1]
    ang2 = md.geometry.compute_phi(traj_corrected, opt=True, periodic=False)[1]
    assert not (np.max(np.abs(ang1 - ang2)) < epsilon)


def test_angle_pbc_fails1(get_fn):
    traj_uncorrected = md.load(get_fn('1am7_uncorrected.xtc'), top=get_fn('1am7_protein.pdb'))
    indices = []
    for i in range(traj_uncorrected.n_residues):
        r = traj_uncorrected.topology.residue(i)
        indices.append((r.atom(0).index, r.atom(1).index, r.atom(2).index))
    traj_corrected = md.load(get_fn('1am7_corrected.xtc'), top=get_fn('1am7_protein.pdb'))

    epsilon = 1E-2

    ang1 = md.geometry.compute_angles(traj_uncorrected, indices, opt=False, periodic=False)
    ang2 = md.geometry.compute_angles(traj_corrected, indices, opt=False, periodic=False)
    assert not (np.max(np.abs(ang1 - ang2)) < epsilon)


def test_no_indices(get_fn):
    for fn in ['2EQQ.pdb', '1bpi.pdb']:
        for opt in [True, False]:
            t = md.load(get_fn(fn))
            assert md.compute_distances(t, np.zeros((0, 2), dtype=int), opt=opt).shape == (t.n_frames, 0)
            assert md.compute_angles(t, np.zeros((0, 3), dtype=int), opt=opt).shape == (t.n_frames, 0)
            assert md.compute_dihedrals(t, np.zeros((0, 4), dtype=int), opt=opt).shape == (t.n_frames, 0)
