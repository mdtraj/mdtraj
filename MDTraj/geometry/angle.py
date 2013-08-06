# This file is part of MDTraj.
#
# Copyright 2013 Stanford University
#
# MSMBuilder is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

##############################################################################
# Imports
##############################################################################

import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import _HAVE_OPT
if _HAVE_OPT:
    from mdtraj.geometry import ffi, C
    from mdtraj.utils.ffi import cpointer

##############################################################################
# Functions
##############################################################################

def compute_angles(traj, angle_indices, opt=True):
    """Compute the bond angles between the supplied triplets of indices
    for each frame in a trajectory

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    angle_indices : np.ndarray, shape=(num_pairs, 2), dtype=int
       Each row gives the indices of three atoms which together make an angle.
    opt : bool, default=True
        Use an optimized native library to calculate distances. Using this
        library requires the python package "cffi" (c foreign function
        interface) which is installable via "easy_install cffi" or "pip
        install cffi". See https://pypi.python.org/pypi/cffi for more details.

    Returns
    -------
    angles : np.ndarray, shape=[n_frames, n_angles], dtype=float
        The angles are in radians
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3))
    triplets = ensure_type(np.asarray(angle_indices), dtype=np.int32, ndim=2, name='angle_indices', shape=(None, 2))
    out = np.zeros((xyz.shape[0], triplets.shape[0]), dtype=np.float32)
    if _HAVE_OPT and opt:
        C.angles(cpointer(xyz), cpointer(triplets), cpointer(out), xyz.shape[0],
                 xyz.shape[1], triplets.shape[0])
    else:
        _angles(xyz, triplets, out)
    return out


def _angles(xyz, angle_indices, out):
    for i in xrange(xyz.shape[1]):
        for j, (m, o, n) in enumerate(angle_indices):
            u_prime = xyz[i, m, :] - xyz[i, o, :]
            v_prime = xyz[i, n, :] - xyz[i, o, :]
            u_norm = np.linalg.norm(u_prime)
            v_norm = np.linalg.norm(v_prime)

            out[i, j] = np.arccos(np.dot(u_prime / u_norm, v_prime / v_norm))

    return angles

