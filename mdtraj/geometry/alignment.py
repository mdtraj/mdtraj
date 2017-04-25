##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Kyle A. Beachamp
# Contributors: Robert McGibbon
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
Here, we have implemented RMSD calculation using both Kabsch algorithm and
via the quaterion-based characteristic polynomial.  Both are implemented
in pure python/numpy.

Notes
-----
This module works at a frame-by-frame level, not the trajectory level.

Reference
---------
Rapid calculation of RMSDs using a quaternion-based characteristic polynomial.
Acta Crystallogr A 61(4):478-480.
"""

from __future__ import print_function, division
import numpy as np
import scipy.optimize
from mdtraj.utils import ensure_type


class Transformation(object):
    """Operator capable of rotating and translating a conformation

    Parameters
    ----------
    rotation : np.array, shape=(3,3)
        Rotation matrix
    translation : np.array, shape=(3)
        Translation vector
    """
    def __init__(self, rotation, translation):
        self.rotation = rotation
        self.translation = translation

    def transform(self, coordinates):
        """Apply an affine transformation to a set of coordinates.

        Parameters
        ----------
        coordinates : ndarray, shape = (n_atoms, 3)
            xyz coordinates of a `single` frame.

        Returns
        -------
        mapped_coordinates : ndarray, shape = (n_atoms, 3)
            xyz coordinates after being rotated and translated

        """
        return coordinates.dot(self.rotation) + self.translation


    def __call__self(self, coordinates):
        return self.transform(coordinates)


def compute_transformation(mobile, target):
    """Returns a Transformation object to align mobile onto target.

    Parameters
    ----------
    mobile : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame, to be aligned onto target.
    target : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame

    Returns
    -------
    T : Transformation aligning mobile to target.
    """
    translation, rotation = compute_translation_and_rotation(mobile, target)
    return Transformation(rotation, translation)


def transform(mobile, target):
    """Align mobile onto target and return transformed coordinates.

    Parameters
    ----------
    mobile : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame, to be aligned onto target.
    target : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame

    Returns
    -------
    mobile_prime : ndarray, shape = (n_atoms, 3)
        Transformed coordinates of mobile, optimally aligned to target.
    """
    T = compute_transformation(mobile, target)
    mobile_prime = T.transform(mobile)
    return mobile_prime


def compute_translation_and_rotation(mobile, target):
    """Returns the translation and rotation mapping mobile onto target.

    Parameters
    ----------
    mobile : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame, to be aligned onto target.
    target : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame

    Returns
    -------
    translation : ndarray, shape=(3,)
        Difference between the centroids of the two conformations
    rotation : ndarray, shape=(3,3)
        Rotation matrix to apply to mobile to carry out the transformation.
    """

    ensure_type(mobile, 'float', 2, 'mobile', warn_on_cast=False, shape=(None, 3))
    ensure_type(target, 'float', 2, 'target', warn_on_cast=False, shape=(target.shape[0], 3))

    mu1 = mobile.mean(0)
    mu2 = target.mean(0)

    mobile = mobile - mu1
    target = target - mu2

    correlation_matrix = np.dot(np.transpose(mobile), target)
    V, S, W_tr = np.linalg.svd(correlation_matrix)
    is_reflection = (np.linalg.det(V) * np.linalg.det(W_tr)) < 0.0
    if is_reflection:
        V[:, -1] = -V[:, -1]
    rotation = np.dot(V, W_tr)

    translation = mu2 - mu1.dot(rotation)

    return translation, rotation


def rmsd_kabsch(xyz1, xyz2):
    """Returns the RMSD distance between conformations xyz1 and xyz2.

    Parameters
    ----------
    xyz1 : ndarray, shape = (n_atoms, 3)
        xyz1 coordinates of a SINGLE frame
    xyz2 : ndarray, shape = (n_atoms, 3)
        xyz2 coordinates of a SINGLE frame

    Returns
    -------
    rmsd : float
        RMSD between xyz1 and xyz2
    """
    xyz1_prime = transform(xyz1, xyz2)
    delta = xyz1_prime - xyz2
    rmsd = (delta ** 2.0).sum(1).mean() ** 0.5
    return rmsd


def _center(conformation):
    """Center the conformation"""
    ensure_type(conformation, 'float', 2, 'conformation', warn_on_cast=False, shape=(None, 3))
    centroid = np.mean(conformation, axis=0)
    centered = conformation - centroid
    return centered


def rmsd_qcp(conformation1, conformation2):
    """Compute the RMSD with Theobald's quaterion-based characteristic
    polynomial

    Rapid calculation of RMSDs using a quaternion-based characteristic polynomial.
    Acta Crystallogr A 61(4):478-480.

    Parameters
    ----------
    conformation1 : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the first conformation
    conformation2 : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the second conformation

    Returns
    -------
    rmsd : float
        The root-mean square deviation after alignment between the two pointsets
    """
    ensure_type(conformation1, np.float32, 2, 'conformation1', warn_on_cast=False, shape=(None, 3))
    ensure_type(conformation2, np.float32, 2, 'conformation2', warn_on_cast=False, shape=(conformation1.shape[0], 3))

    A = _center(conformation1)
    B = _center(conformation2)
    if not A.shape[0] == B.shape[0]:
        raise ValueError('conformation1 and conformation2 must have same number of atoms')
    n_atoms = len(A)

    # the inner product of the structures A and B
    G_A = np.einsum('ij,ij', A, A)
    G_B = np.einsum('ij,ij', B, B)
    # print 'GA', G_A, np.trace(np.dot(A.T, A))
    # print 'GB', G_B, np.trace(np.dot(B.T, B))

    # M is the inner product of the matrices A and B
    M = np.dot(B.T, A)

    # unpack the elements
    Sxx, Sxy, Sxz = M[0, :]
    Syx, Syy, Syz = M[1, :]
    Szx, Szy, Szz = M[2, :]

    # do some intermediate computations to assemble the characteristic
    # polynomial
    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    # two of the coefficients
    C2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C1 = 8.0 * (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx - Sxx * Syy * Szz - Syz * Szx * Sxy - Szy * Syx * Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    # the other coefficient
    C0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 \
        + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) \
        + (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz)) * (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz)) \
        + (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz)) * (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz)) \
        + (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz)) * (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz)) \
        + (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz)) * (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz))

    E0 = (G_A + G_B) / 2.0
    f = lambda x: x ** 4.0 + C2 * x ** 2. + C1 * x + C0
    df = lambda x: 4 * x ** 3.0 + 2 * C2 * x + C1
    max_eigenvalue = scipy.optimize.newton(f, E0, df)
    rmsd = np.sqrt(np.abs(2.0 * (E0 - max_eigenvalue) / n_atoms))
    return rmsd

def compute_average_structure(xyz):
    """Compute the average structure from a set of frames.

    The frames are first aligned to minimize the RMSD from the average structure,
    and then the average is returned.

    Parameters
    ----------
    xyz : ndarray, shape = (n_frames, n_atoms, 3)
        xyz coordinates of each atom in each frame

    Returns
    -------
    average : ndarray, shape = (n_atoms, 3)
        the average structure
    """
    n_frames = xyz.shape[0]
    n_atoms = xyz.shape[1]
    candidate = xyz[0]

    # In practice there is usually negligible improvement after the second iteration.
    # We do three just to be safe.

    for iteration in range(3):
        average = np.zeros((n_atoms, 3))
        for frame in range(n_frames):
            average += transform(xyz[frame], candidate)
        candidate = average/n_frames
    return candidate
