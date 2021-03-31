##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2016 Stanford University and the Authors
#
# Authors: Peter Eastman
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


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import numpy as np
from mdtraj.utils import ensure_type

from libcpp.vector cimport vector

__all__ = ['compute_neighborlist']

cdef extern from "neighborlist.h":
    vector[vector[int]] _compute_neighborlist(float* positions, int n_atoms, float cutoff, float* box_matrix) nogil


##############################################################################
# Functions
##############################################################################

def compute_neighborlist(traj, cutoff, frame=0, periodic=True):
    """compute_neighborlist(traj, cutoff, frame=0, periodic=True)

    Find (spatially) neighboring atoms in a trajectory.

    For every atom, build a list of all other atoms that are within a cutoff
    distance of it.

    Parameters
    ----------
    traj : md.Trajectory
        An MDTraj trajectory
    cutoff : float
        Distance cutoff to define 'neighboring'
    frame : int
        Index of the frame of the trajectory to build the neighbor list based on
    periodic : bool
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.

    Returns
    -------
    neighbors : list of np.ndarray, dtype=int
        List of arrays, of length n_atoms. neighbors[i] contains the indices of all
        atoms within the cutoff distance of atom i.
    """

    cdef float[:, ::1] xyz = traj.xyz[frame]
    cdef float[:, ::1] box_matrix
    cdef float* box_matrix_pointer
    cdef int is_periodic = periodic and (traj.unitcell_vectors is not None)
    if is_periodic:
        unitcell_vectors = ensure_type(traj.unitcell_vectors[frame],
                                       dtype=np.float32, ndim=2,
                                       name='unitcell_vectors',
                                       warn_on_cast=False)
        box_matrix = np.asarray(unitcell_vectors, order='c')
        box_matrix_pointer = &box_matrix[0,0]
    else:
        box_matrix_pointer = NULL
    cdef int n_atoms = traj.xyz.shape[1]
    cdef vector[vector[int]] atom_neighbors = _compute_neighborlist(&xyz[0,0], n_atoms, cutoff, box_matrix_pointer)

    # Copy the results over to Numpy arrays.

    cdef int[::1] atom_neighbors_mview
    cdef list neighbors = []
    cdef int i
    for i in range(n_atoms):
        if atom_neighbors[i].size() > 0:
            atom_neighbors_mview = <int[:atom_neighbors[i].size()]> (<int*> (&atom_neighbors[i][0]))
            neighbors.append(np.array(atom_neighbors_mview, dtype=int, copy=True))
        else:
            neighbors.append(np.empty(0, dtype=int))
    return neighbors
