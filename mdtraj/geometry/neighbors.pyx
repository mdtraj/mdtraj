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

cdef extern from "neighbors.hpp":
    vector[int] _compute_neighbors(float* xyz, int n_atoms, float cutoff,
        vector[int]& query_indices, vector[int]& haystack_indices,
        float* box_matrix) nogil

##############################################################################
# Functions
##############################################################################

def compute_neighbors(traj, cutoff, query_indices, haystack_indices=None,
                      periodic=True):
    """compute_neighbors(traj, cutoff, query_indices, haystack_indices=None, periodic=True)

    Find (spatially) neighboring atoms in a trajectory.

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
    haystack_indices : np.ndarray, shape=(n_query_indices,), dtype=int, optional
        If supplied, restrict the search to only those atoms in
        `haystack_indices`.
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

    query_indices = ensure_type(query_indices, dtype=np.int32, ndim=1,
                                name='query_indices', warn_on_cast=False)
    haystack_indices = ensure_type(haystack_indices, dtype=np.int32, ndim=1,
                                   name='haystack_indices', warn_on_cast=False,
                                   can_be_none=True)
    if haystack_indices is None:
        haystack_indices = np.arange(traj.xyz.shape[1])

    if not np.all((query_indices >= 0) * (query_indices < traj.xyz.shape[1]) * (query_indices < traj.xyz.shape[1])):
        raise ValueError("query_indices must be valid positive indices")
    if not np.all((haystack_indices >= 0) * (haystack_indices < traj.xyz.shape[1]) * (haystack_indices < traj.xyz.shape[1])):
        raise ValueError("haystack_indices must be valid positive indices")

    cdef int i
    cdef int n_frames = traj.xyz.shape[0]
    cdef float[:, :, ::1] xyz = traj.xyz
    cdef float[:, :, ::1] box_matrix
    cdef float* box_matrix_pointer
    cdef vector[int] query_indices_ = query_indices
    cdef vector[int] haystack_indices_ = haystack_indices
    cdef vector[int] frame_neighbors
    cdef int[::1] frame_neighbors_mview
    cdef int is_periodic = periodic and (traj.unitcell_vectors is not None)
    if is_periodic:
        unitcell_vectors = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', warn_on_cast=False)
        box_matrix = np.asarray(unitcell_vectors, order='c')
    else:
        box_matrix_pointer = NULL

    results = []  # list of numpy arrays
    for i in range(n_frames):
        if is_periodic:
            box_matrix_pointer = &box_matrix[i,0,0]
        frame_neighbors = _compute_neighbors(
            &xyz[i,0,0], traj.xyz.shape[1], cutoff, query_indices_,
            haystack_indices_, box_matrix_pointer)
        # now, we need to go from STL vector[int] to a numpy array without
        # egregious copying performance.
        # I can't find any great cython docs on this...
        # first, convert to a memoryview by casting the pointer to the memory
        # block. this implements the python buffer protocol...
        if frame_neighbors.size() > 0:
            frame_neighbors_mview = <int[:frame_neighbors.size()]> (<int*> (&frame_neighbors[0]))
        # so then we can copy it into numpy-managed memory. copy is necessary,
        # because once we exit this scope, C++ will clean up the vector.
            results.append(np.array(frame_neighbors_mview, dtype=int, copy=True))
        else:
            results.append(np.empty(0, dtype=int))

    return results

