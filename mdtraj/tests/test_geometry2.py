import numpy as np
import scipy.spatial.distance

import mdtraj as md
from mdtraj.geometry import geometry
from mdtraj.testing import get_fn
from mdtraj.testing.cpptraj import cpptraj_distance

def test_3():
    t = md.load(get_fn('tip3p_300K_1ATM.pdb'))
    t = t.atom_slice([0, 1])

    t.unitcell_lengths = np.array([[2., 2., 2.]])
    # t.unitcell_angles = np.array([[ 77.84587097,  22.20738029,  88.23445892]])
    t.unitcell_angles = np.array([[ 16.22463799,  28.89242363,  45.09446335]])
    t.xyz = np.array([[[ 1.89890003,  1.05920005,  0.33930001],
                       [ 1.94840002,  1.21379995,  1.35549998]]])

    print('\ncpptraj     ', cpptraj_distance(t, 0, 1))
    print('mdtraj1     ', md.compute_distances(t, [[0, 1]], opt=False))
    print('brute force ', simple_brute_force(t.xyz[0,0], t.xyz[0,1],
          t.unitcell_vectors[0,0],
          t.unitcell_vectors[0,1],
          t.unitcell_vectors[0,2]))
    print('mdtraj2     ', geometry.compute_distances(
               t.xyz,
               np.asarray(t.unitcell_vectors, order='c'),
               np.array([[0,1]])))
    # print('openmm      ', openmm_distance(t))


def test_2():
    t = md.load(get_fn('tip3p_300K_1ATM.pdb'))
    t = t.atom_slice([0, 1])

    t.unitcell_lengths = np.array([[2., 2., 2.]])
    t.unitcell_angles = np.array([[ 46.77410889,  63.44205856,  47.81904221]])
    t.xyz = np.array([[[ 0.85960001,  1.70850003,  1.29089999], [ 0.8581,      1.45210004,  0.59310001]]])

    print('\ncpptraj     ', cpptraj_distance(t, 0, 1))
    print('mdtraj1     ', md.compute_distances(t, [[0, 1]], opt=False))
    print('brute force ', simple_brute_force(t.xyz[0,0], t.xyz[0,1],
          t.unitcell_vectors[0,0],
          t.unitcell_vectors[0,1],
          t.unitcell_vectors[0,2]))
    print('mdtraj2     ', geometry.compute_distances(
               t.xyz,
               np.asarray(t.unitcell_vectors, order='c'),
               np.array([[0,1]])))
    # print('openmm      ', openmm_distance(t))


def simple_brute_force(v1, v2, a, b, c):
    c = c - b*np.round(c[1]/b[1])
    c = c - a*np.round(c[0]/a[0])
    b = b - a*np.round(b[0]/a[0])

    H = np.array([a,b,c]).T
    Hi = np.linalg.inv(H)

    f1 = np.dot(Hi, v1)
    f2 = np.dot(Hi, v2)

    f1 = f1 - np.floor(f1)
    f2 = f2 - np.floor(f2)

    v1 = np.dot(H, f1)
    v2 = np.dot(H, f2)

    d = float('infinity')
    for i in range(-5, 5):
        for j in range(-5, 5):
            for k in range(-5, 5):
                v2_p = v2 + (i*a + j*b + k*c)
                d_p = scipy.spatial.distance.euclidean(v1, v2_p)
                if d_p < d:
                    d = d_p
                    qq = (i,j,k)

    return d


def openmm_distance(t, atom1=0, atom2=1):
    t = t.atom_slice([atom1, atom2])
    import simtk.openmm as mm
    from simtk.openmm.app.internal.unitcell import computePeriodicBoxVectors

    system = mm.System()
    force = mm.CustomNonbondedForce('r')
    force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    for _ in range(2):
        system.addParticle(1)
        force.addParticle([])


    system.addForce(force)
    system.setDefaultPeriodicBoxVectors(*computePeriodicBoxVectors(
        t.unitcell_lengths[0,0],
        t.unitcell_lengths[0,1],
        t.unitcell_lengths[0,2],
        np.deg2rad(t.unitcell_angles[0,0]),
        np.deg2rad(t.unitcell_angles[0,1]),
        np.deg2rad(t.unitcell_angles[0,2]),
    ))
    a,b,c = system.getDefaultPeriodicBoxVectors()
    force.setCutoffDistance(min(a[0]/2, b[1]/2, c[2]/2))

    context = mm.Context(system, mm.VerletIntegrator(0.002), mm.Platform.getPlatformByName('Reference'))
    xyz = t.xyz[0]

    assert system.usesPeriodicBoundaryConditions()

    # a,b,c = system.getDefaultPeriodicBoxVectors()
    # xyz[1] = xyz[1] + 100*a._value

    context.setPositions(xyz)
    return context.getState(getEnergy=True).getPotentialEnergy()._value


def test_1():
    t = md.load(get_fn('tip3p_300K_1ATM.pdb'))

    random = np.random.RandomState(0)

    for i in range(500):
        t.unitcell_angles = random.uniform(low=15.0, high=90.0, size=(1,3))
        if not np.all(np.isfinite(t.unitcell_vectors)):
            continue

        pairs = random.randint(low=0, high=t.n_atoms-1, size=(1,2))

        print('angles', repr(t.unitcell_angles[0]))
        print('lengths', repr(t.unitcell_lengths[0]))
        print('coors', repr(t.xyz[0,pairs[0,0]]), repr(t.xyz[0,pairs[0,1]]))

        d1 = geometry.compute_distances(
            t.xyz,
            np.asarray(t.unitcell_vectors, order='c'),
            pairs)[0,0]

        d2 = simple_brute_force(v1=t.xyz[0, pairs[0,0]],
                         v2=t.xyz[0, pairs[0,1]],
                         a=t.unitcell_vectors[0, 0],
                         b=t.unitcell_vectors[0, 1],
                         c=t.unitcell_vectors[0, 2])

        d3 = cpptraj_distance(t, pairs[0,0], pairs[0,1])[0]
        # d4 = openmm_distance(t, pairs[0,0], pairs[0,1])
        print('new libgeometry', d1)
        print('brute force    ', d2)
        print('cpptraj        ', d3)
        # print('openmm         ', d4)

        # np.testing.assert_almost_equal(d1, d3, decimal=4)
        np.testing.assert_almost_equal(d1, d2, decimal=4)
        # np.testing.assert_almost_equal(d1, d4, decimal=4)


