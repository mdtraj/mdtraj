import numpy as np
from mdtraj.testing import eq
from mdtraj.utils import uniform_quaternion, rotation_matrix_from_quaternion

def test_quaternion_0():
    q = uniform_quaternion(random_state=2)
    eq(np.sum(q**2), np.float64(1.0))
    eq(q.shape, (4,))

    q2 = uniform_quaternion(size=2)
    eq(np.sum(q2**2, axis=1), np.ones(2))
    eq(q2.shape, (2,4))

    q2 = uniform_quaternion(size=(6,6))
    eq(np.sum(q2**2, axis=2), np.ones((6,6)))
    eq(q2.shape, (6,6,4))


def test_quaternion_1():
    q = uniform_quaternion(random_state=2)
    rot = rotation_matrix_from_quaternion(q)
    
    eq(rot.shape, (3,3))
    eq(np.linalg.det(rot), np.double(1.0))
    eq(np.linalg.inv(rot), rot.T)


def test_quaternion_2():
    q = uniform_quaternion(size=2)
    rot = rotation_matrix_from_quaternion(q)
    
    eq(rot.shape, (2, 3, 3))

    # check 1st rotation matrix
    eq(np.linalg.det(rot[0]), np.double(1.0))
    eq(np.linalg.inv(rot[0]), rot[0].T)
    # check 2nd rotation matrix
    eq(np.linalg.det(rot[1]), np.double(1.0))
    eq(np.linalg.inv(rot[1]), rot[1].T)
