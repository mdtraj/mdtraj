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

from libcpp.vector cimport vector

__all__ = ['compute_neighbors']

cdef extern:
    vector[int] _compute_neighbors(float* xyz, int n_atoms, float cutoff,
        vector[int]& query_indices, vector[int]& haystack_indices,
        float* box_matrix) nogil

##############################################################################
# Functions
##############################################################################

def compute_neighbors(traj, cutoff, query_indices, haystack_indices=None,
                      periodic=True):
    """Find (spatially) neighboring atoms in a trajectory.

    Given a set of query_indices representing and a distance cutoff, compute
    the indices of all atoms whose distance to 1 or more of the query points
    is less than cutoff.

    Parameters
    ----------
    traj : md.Trajectory
        An MDTraj trajectory
    cutoff : float
        Distance cutoff to define 'neighboring'
    query_indices : np.ndarray, shape=(n_query_indices,), dtype=int
        The matching atoms are those that are within `cutoff` of one or more
        of the atoms with indices in `query_indices`.
    haystack_indies : np.ndarray, shape=(n_query_indices,), dtype=int, optional
        If supplied, restrict the search to only those atoms in
        `haystack_indies`.
    periodic : bool
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.

    Returns
    -------
    matches : list of np.ndarray, shape=(n_matches,), dtype=int
        List of arrays, of length n_frames. Each item in the list is a 1D array
        of the indices of the matching atoms.
    """
    if traj.topology is None:
        raise ValueError('traj must have a topology defined')
    else:
        query_indices = ensure_type(query_indices, dtype=np.int32, ndim=1, name='query_indices', warn_on_cast=False)
        haystack_indices = ensure_type(haystack_indices, dtype=np.int32, ndim=1, name='haystack_indices', warn_on_cast=False)
        if haystack_indices is None:
            haystack_indices = np.arange(traj.xyz.shape[1])

        if not np.all((query_indices >= 0) * (query_indices < traj.xyz.shape[1]) * (query_indices < traj.xyz.shape[1])):
            raise ValueError("query_indices must be valid positive indices")
        if not np.all((haystack_indices >= 0) * (haystack_indices < traj.xyz.shape[1]) * (haystack_indices < traj.xyz.shape[1])):
            raise ValueError("haystack_indices must be valid positive indices")

    cdef float[:, :, ::1] xyz = traj.xyz
    cdef float[:, :, ::1] box_matrix
    if traj.unitcell_vectors is not None:
        box_matrix = traj.unitcell_vectors

    return _compute_neighbors(
        &xyz[0,0,0], traj.xyz.n_atoms, cutoff, query_indices,
        haystack_indices, &box_matrix[0,0,0])
