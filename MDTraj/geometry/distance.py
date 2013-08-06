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
try:
    from mdtraj.geometry._distance import distance as _opt_distance
    from mdtraj.geometry._distance import displacement as _opt_displacement
    _HAVE_OPT = True
except ImportError:
    _HAVE_OPT = False


##############################################################################
# Functions
##############################################################################


def compute_distances(traj, atom_pairs, periodic=True):
    """Compute the distances between pairs of atoms in each frame.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    atom_pairs : np.ndarray, shape=(num_pairs, 2), dtype=int
        Each row gives the indices of two atoms involved in the interaction.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.

    Returns
    -------
    distances : np.ndarray, shape=(n_frames, num_pairs), dtype=float
        The distance, in each frame, between each pair of atoms.
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='taj.xyz', shape=(None, None, 3))
    pairs = ensure_type(np.asarray(atom_pairs), dtype=np.int32, ndim=2, name='atom_pairs', shape=(None, 2))

    if periodic is True and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(xyz), 3, 3))
        if _HAVE_OPT:
            return _opt_distance(xyz, pairs, box)
        return _distance_mic(xyz, pairs, box)

    # either there are no unitcell vectors or they dont want to use them
    if _HAVE_OPT:
        return _opt_distance(xyz, pairs)
    return _distance_mic(yz, pairs)


def compute_displacements(traj, atom_pairs, periodic=True):
    """Compute the displacement vector between pairs of atoms in each frame

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute distances in
    atom_pairs : np.ndarray, shape[num_pairs, 2], dtype=int
        Each row gives the indices of two atoms.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.

    Returns
    -------
    displacements : np.ndarray, shape=[n_frames, n_pairs, 3], dtype=float32
         The displacememt vector, in each frame, between each pair of atoms.
    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='taj.xyz', shape=(None, None, 3))
    pairs = ensure_type(np.asarray(atom_pairs), dtype=np.int32, ndim=2, name='atom_pairs', shape=(None, 2))

    if periodic is True and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(xyz), 3, 3))
        if _HAVE_OPT:
            return _opt_displacement(xyz, pairs, box)
        return _displacement_mic(xyz, pairs, box)

    # either there are no unitcell vectors or they dont want to use them
    if _HAVE_OPT:
        return _opt_displacement(xyz, pairs)
    return _displacement(xyz, pairs)


##############################################################################
# pure python implementation of the core routines
##############################################################################

def _distance(xyz, pairs):
    "Distance between pairs of points in each frame"
    delta = np.diff(xyz[:, pairs], axis=2)[:, :, 0]
    return (delta ** 2.).sum(-1) ** 0.5


def _displacement(xyz, pairs):
    "Displacement vector between pairs of points in each frame"
    return np.diff(xyz[:, pairs], axis=2)[:, :, 0]


def _distance_mic(xyz, pairs, box_vectors):
    """Distance between pairs of points in each frame under the minimum image
    convention for periodic boundary conditions.

    The computation follows scheme B.9 in Tukerman, M. "Statistical
    Mechanics: Theory and Molecular Simulation", 2010.

    This is a slow pure python implementation, mostly for testing.
    """
    out = np.empty((xyz.shape[0], pairs.shape[0]), dtype=np.float32)
    for i in range(len(xyz)):
        hinv = np.linalg.inv(box_vectors[i])

        for j, (a,b) in enumerate(pairs):
            s1 = np.dot(hinv, xyz[i,a,:])
            s2 = np.dot(hinv, xyz[i,b,:])
            s12 = s2 - s1

            s12 = s12 - np.round(s12)
            r12 = np.dot(box_vectors[i], s12)
            out[i, j] = np.sqrt(np.sum(r12 * r12))
    return out


def _displacement_mic(xyz, pairs, box_vectors):
    """Displacement vector between pairs of points in each frame under the
    minimum image convention for periodic boundary conditions.

    The computation follows scheme B.9 in Tukerman, M. "Statistical
    Mechanics: Theory and Molecular Simulation", 2010.

    This is a slow pure python implementation, mostly for testing.
    """
    out = np.empty((xyz.shape[0], pairs.shape[0], 3), dtype=np.float32)
    for i in range(len(xyz)):
        hinv = np.linalg.inv(box_vectors[i])

        for j, (a,b) in enumerate(pairs):
            s1 = np.dot(hinv, xyz[i,a,:])
            s2 = np.dot(hinv, xyz[i,b,:])
            s12 = s2 - s1
            s12 = s12 - np.round(s12)
            out[i, j] = np.dot(box_vectors[i], s12)

    return out
