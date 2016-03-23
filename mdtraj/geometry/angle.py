##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A Beauchamp
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

from __future__ import print_function, division
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import _geometry, distance
import warnings

__all__ = ['compute_angles']

##############################################################################
# Functions
##############################################################################


def compute_angles(traj, angle_indices, periodic=True, opt=True):
    """Compute the bond angles between the supplied triplets of indices in each frame of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mdtraj trajectory.
    angle_indices : np.ndarray, shape=(num_angles, 3), dtype=int
       Each row gives the indices of three atoms which together make an angle.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat angles that cross periodic images using
        the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate distances. Our optimized
        SSE angle calculation implementation is 10-20x faster than the
        (itself optimized) numpy implementation.

    Returns
    -------
    angles : np.ndarray, shape=[n_frames, n_angles], dtype=float
        The angles are in radians
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    triplets = ensure_type(angle_indices, dtype=np.int32, ndim=2, name='angle_indices', shape=(None, 3), warn_on_cast=False)
    if not np.all(np.logical_and(triplets < traj.n_atoms, triplets >= 0)):
        raise ValueError('angle_indices must be between 0 and %d' % traj.n_atoms)

    if len(triplets) == 0:
        return np.zeros((len(xyz), 0), dtype=np.float32)

    out = np.zeros((xyz.shape[0], triplets.shape[0]), dtype=np.float32)
    if periodic is True and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(xyz), 3, 3))
        if opt:
            orthogonal = np.allclose(traj.unitcell_angles, 90)
            _geometry._angle_mic(xyz, triplets, box.transpose(0, 2, 1).copy(), out, orthogonal)
            return out
        else:
            _angle(traj, triplets, periodic, out)
            return out

    if opt:
        _geometry._angle(xyz, triplets, out)
    else:
        _angle(traj, triplets, periodic, out)
    return out


def _angle(traj, angle_indices, periodic, out):
    ix01 = angle_indices[:, [1, 0]]
    ix21 = angle_indices[:, [1, 2]]

    u_prime = distance.compute_displacements(traj, ix01, periodic=periodic, opt=False)
    v_prime = distance.compute_displacements(traj, ix21, periodic=periodic, opt=False)
    u_norm = np.sqrt((u_prime**2).sum(-1))
    v_norm = np.sqrt((v_prime**2).sum(-1))

    # adding a new axis makes sure that broasting rules kick in on the third
    # dimension
    u = u_prime / (u_norm[..., np.newaxis])
    v = v_prime / (v_norm[..., np.newaxis])

    return np.arccos((u * v).sum(-1), out=out)
