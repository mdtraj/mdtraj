##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Christoph Klein
# Contributors: Tim Moore
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

from __future__ import print_function, division

import numpy as np

from mdtraj.geometry.distance import compute_center_of_geometry

__all__ = ['compute_gyration_tensor', 'principal_moments', 'asphericity', 'acylindricity', 'relative_shape_antisotropy', 'relative_shape_anisotropy']


def compute_gyration_tensor(traj):
    """Compute the gyration tensor of a trajectory.

    For every frame,

    .. math::

        S_{xy} = sum_{i_atoms} r^{i}_x r^{i}_y

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute gyration tensor of.

    Returns
    -------
    S_xy:  np.ndarray, shape=(traj.n_frames, 3, 3), dtype=float64
        Gyration tensors for each frame.
    
    References
    ----------
    .. [1] https://isg.nist.gov/deepzoomweb/measurement3Ddata_help#shape-metrics-formulas

    """
    center_of_geom = np.expand_dims(compute_center_of_geometry(traj), axis=1)
    xyz = traj.xyz - center_of_geom
    return np.einsum('...ji,...jk->...ik', xyz, xyz) / traj.n_atoms

def principal_moments(traj):
    """Compute the principal moments of a trajectory.

    For each frame calculate the eigenvalues of the gyration tensor.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute gyration tensor of.

    Returns
    -------
    evecs:  np.ndarray, shape=(traj.n_frames, 3), dtype=float64
        Principal moments for each frame in ascending order.

    """
    gyration_tensor = compute_gyration_tensor(traj)
    return np.linalg.eigvalsh(gyration_tensor)

def asphericity(traj):
    """Compute the asphericity of a trajectory.

    For each frame compute the principal moments then,

    .. math::

        b = \frac{1}{2}(\lambda_1^2 + \lambda_2^2)


    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute gyration tensor of.

    Returns
    -------
    b:  np.ndarray, shape=(traj.n_frames, 1), dtype=float64
        Asphericity of each frame of the trajectory.

    """
    pm = principal_moments(traj)
    b = pm[:,2] - (pm[:,0] + pm[:,1]) / 2.0
    return b

def acylindricity(traj):
    """Compute the acylindricity of a trajectory.

    For each frame compute the principal moments then,

    .. math::

        c = \lambda_2^2 - \lambda_1^2


    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute gyration tensor of.

    Returns
    -------
    c:  np.ndarray, shape=(traj.n_frames, 1), dtype=float64
        Acylindricity of each frame of the trajectory.

    """
    pm = principal_moments(traj)
    c = pm[:,1] - pm[:,0]
    return c

def relative_shape_anisotropy(traj):
    """Compute the relative shape anisotropy of a trajectory.

    For each frame compute the principal moments then,

    .. math::

        \kappa^2 = \frac{3}{2}\frac{\lambda_1^4 + \lambda_2^4 + \lambda_3^4}{(\lambda_1^2 + \lambda_2^2 + \lambda_3^2)^2} - \frac{1}{2}


    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute gyration tensor of.

    Returns
    -------
    c:  np.ndarray, shape=(traj.n_frames, 1), dtype=float64
        Relative shape anisotropy of each frame of the trajectory.

    """
    pm = principal_moments(traj)
    kappa2 = 1.5 * np.square(pm).sum(axis=1) / np.square(pm.sum(axis=1)) - 0.5
    return kappa2

relative_shape_antisotropy = relative_shape_anisotropy

def _compute_gyration_tensor_slow(traj):
    """Compute the gyration tensor of a trajectory. """
    xyz = traj.xyz
    center_of_geom = np.expand_dims(compute_center_of_geometry(traj), axis=1)
    centered_xyz = xyz - center_of_geom

    S_nm = np.zeros(shape=(traj.n_frames, 3, 3), dtype=np.float64)
    for n, xyz in enumerate(centered_xyz):
        N = xyz.shape[0]
        for r in xyz:
            S_nm[n, 0, 0] += r[0] * r[0]
            S_nm[n, 1, 1] += r[1] * r[1]
            S_nm[n, 2, 2] += r[2] * r[2]
            S_nm[n, 0, 1] += r[0] * r[1]
            S_nm[n, 0, 2] += r[0] * r[2]
            S_nm[n, 1, 2] += r[1] * r[2]
        S_nm[n, 1, 0] = S_nm[n, 0, 1]
        S_nm[n, 2, 0] = S_nm[n, 0, 2]
        S_nm[n, 2, 1] = S_nm[n, 1, 2]
        S_nm[n, :, :] /= N
    return S_nm