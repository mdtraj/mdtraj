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

__all__ = ['compute_gyration_tensor']


def compute_gyration_tensor(traj):
    """Compute the gyration tensor of a trajectory.

    For each frame,

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