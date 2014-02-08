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
from mdtraj.geometry import _geometry

__all__ = ['compute_angles']

##############################################################################
# Functions
##############################################################################


def compute_angles(traj, angle_indices, opt=True):
    """Compute the bond angles between the supplied triplets of indices in each frame of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    angle_indices : np.ndarray, shape=(num_pairs, 2), dtype=int
       Each row gives the indices of three atoms which together make an angle.
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
    triplets = ensure_type(np.asarray(angle_indices), dtype=np.int32, ndim=2, name='angle_indices', shape=(None, 3), warn_on_cast=False)
    if not np.all(np.logical_and(triplets < traj.n_atoms, triplets >= 0)):
        raise ValueError('angle_indices must be between 0 and %d' % traj.n_atoms)

    out = np.zeros((xyz.shape[0], triplets.shape[0]), dtype=np.float32)
    if opt:
        _geometry._angle(xyz, triplets, out)
    else:
        _angle(xyz, triplets, out)
    return out


def _angle(xyz, angle_indices, out):
    #for j, (m, o, n) in enumerate(angle_indices):
    u_prime = xyz[:, angle_indices[:, 0], :] - xyz[:, angle_indices[:, 1], :]
    v_prime = xyz[:, angle_indices[:, 2], :] - xyz[:, angle_indices[:, 1], :]
    u_norm = np.sqrt((u_prime**2).sum(-1))
    v_norm = np.sqrt((v_prime**2).sum(-1))

    # adding a new axis makes sure that broasting rules kick in on the third
    # dimension
    u = u_prime / (u_norm[..., np.newaxis])
    v = v_prime / (v_norm[..., np.newaxis])

    out = np.arccos((u * v).sum(-1), out=out)
