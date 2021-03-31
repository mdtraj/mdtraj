from __future__ import print_function
import numpy as np

import mdtraj as md
from mdtraj.testing import eq
from mdtraj import Trajectory, lprmsd
from mdtraj._lprmsd import _munkres
from mdtraj.utils import rotation_matrix_from_quaternion, uniform_quaternion

random = np.random.RandomState(0)


def test_munkres_0():
    result = _munkres(np.array([[7, 4, 3], [6, 8, 5], [9, 4, 4]], dtype=np.double))
    true = np.array([[0, 0, 1],
                     [1, 0, 0],
                     [0, 1, 0]], dtype=np.int32)
    eq(result, true)


def test_lprmsd_null():
    ref = random.randn(1, 10, 3).astype(np.float32)
    new = np.copy(ref)

    value = lprmsd(Trajectory(xyz=new, topology=None), Trajectory(xyz=ref, topology=None))
    eq(value, np.array([0.0], dtype=np.float32), decimal=3)


def test_lprmsd_0():
    # remap a permutation of all the atoms with no rotation
    ref = random.randn(1, 10, 3).astype(np.float32)
    mapping = random.permutation(10)

    new = ref[:, mapping]

    value = lprmsd(Trajectory(xyz=new, topology=None), Trajectory(xyz=ref, topology=None))
    eq(value, np.array([0.0], dtype=np.float32), decimal=3)


def test_lprmsd_1():
    # resolve a random rotation with no permutation
    ref = random.randn(1, 50, 3).astype(np.float32)
    mapping = np.arange(50)
    rot = rotation_matrix_from_quaternion(uniform_quaternion())
    new = ref[:, mapping].dot(rot)

    value = lprmsd(Trajectory(xyz=new, topology=None), Trajectory(xyz=ref, topology=None), permute_groups=[[]])
    assert value[0] < 1e-2


def test_lprmsd_2():
    # resolve a random rotation with some permutation
    ref = random.randn(1, 50, 3).astype(np.float32)
    # first half of the atoms can permute, last 10 are fixed permutation
    mapping = np.concatenate((random.permutation(10), 10 + np.arange(40)))
    rot = rotation_matrix_from_quaternion(uniform_quaternion())
    new = ref[:, mapping].dot(rot)

    value = lprmsd(Trajectory(xyz=new, topology=None), Trajectory(xyz=ref, topology=None),
                   permute_groups=[np.arange(10)])
    assert value[0] < 1e-2


def test_lprmsd_3(get_fn):
    # resolve rotation and permutation togetehr
    t1 = md.load(get_fn('alanine-dipeptide-explicit.pdb'))[0]
    t2 = md.load(get_fn('alanine-dipeptide-explicit.pdb'))[0]

    h2o_o_indices = [a.index for a in t1.topology.atoms if a.residue.name == 'HOH' and a.element.symbol == 'O'][0:20]
    h2o_h_indices = [a.index for a in t1.topology.atoms if a.residue.name == 'HOH' and a.element.symbol == 'H'][0:20]
    backbone_indices = [a.index for a in t1.topology.atoms if a.element.symbol in ['C', 'N']][:5]

    # scramble two groups of indices
    t2.xyz[:, random.permutation(h2o_o_indices)] = t2.xyz[:, h2o_o_indices]
    t2.xyz[:, random.permutation(h2o_h_indices)] = t2.xyz[:, h2o_h_indices]

    # offset the backbone indices slightly
    t2.xyz[:, backbone_indices] += 0.001 * random.randn(len(backbone_indices), 3)

    # rotate everything
    rot = rotation_matrix_from_quaternion(uniform_quaternion())
    t2.xyz[0].dot(rot)

    print('raw distinguishable indices', backbone_indices)
    atom_indices = np.concatenate((h2o_o_indices, backbone_indices))
    value = md.lprmsd(t2, t1, atom_indices=atom_indices, permute_groups=[h2o_o_indices])
    t1.xyz[:, h2o_o_indices] += random.randn(len(h2o_o_indices), 3)

    print('final value', value)
    assert value[0] < 1e-2


def test_lprmsd_4(get_fn):
    t1 = md.load(get_fn('1bpi.pdb'))
    t1.xyz += 0.05 * random.randn(t1.n_frames, t1.n_atoms, 3)
    t2 = md.load(get_fn('1bpi.pdb'))
    # some random indices
    indices = random.permutation(t1.n_atoms)[:t1.n_atoms - 5]

    got = md.lprmsd(t2, t1, atom_indices=indices, permute_groups=[[]])
    ref = md.rmsd(t2, t1, atom_indices=indices)

    eq(got, ref, decimal=3)


def test_lprmsd_5(get_fn):
    t = md.load(get_fn('frame0.h5'))
    t1 = md.load(get_fn('frame0.h5'))

    r = md.rmsd(t, t1, 0)
    a = md.lprmsd(t, t1, 0, permute_groups=[[]], superpose=True)
    eq(a, r, decimal=3)
