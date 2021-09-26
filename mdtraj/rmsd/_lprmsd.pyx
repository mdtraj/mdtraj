##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Lee-Ping Wang
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
import cython
import numpy as np
from mdtraj.utils import ensure_type
import scipy.spatial.distance

cimport numpy as np
from libcpp.vector cimport vector
from libc.stdio cimport printf
from cpython cimport bool
from cython.parallel cimport prange
np.import_array()
assert sizeof(np.int32_t) == sizeof(int)


##############################################################################
# External c/cpp function declarations
##############################################################################
cdef extern from "include/fancy_index.hpp":
    cdef extern void fancy_index2d(const float* A, int nx, int ny,
        const int* indx, int nindx, const int* indy, int nindy, float* out) nogil
cdef extern float msd_atom_major(int nrealatoms, int npaddedatoms,  float* a,
    float* b, float G_a, float G_b, int computeRot, float rot[9]) nogil
cdef extern float rot_atom_major(const int n_atoms, float* a, const float rot[9]) nogil
cdef extern void sgemm33(const float A[9], const float B[9], float out[9]) nogil
cdef extern void inplace_center_and_trace_atom_major(float* coords, float* traces,
    const int n_frames, const int n_atoms) nogil
cdef extern from "include/Munkres.h":
    cdef cppclass Munkres:
        Munkres()
        void solve(double* icost, int* answer, int m, int n)
cdef extern from "include/euclidean_permutation.hpp":
   vector[int] euclidean_permutation(float* target, float* reference,
       int n_atoms, int n_dims, vector[vector[int]]& permute_groups) nogil
cdef extern from "math.h":
    float sqrtf(float x) nogil


##############################################################################
# Functions
##############################################################################

def lprmsd(target, reference, int frame=0, atom_indices=None, permute_groups=None,
           bool parallel=True, bool superpose=False):
    """lprmsd(target, reference, frame=0, atom_indices=None, permute_groups=None, bool parallel=True, bool superpose=False)

    Compute Linear-Programming Root-Mean-Squared Deviation (LP-RMSD) of all
    conformations in target to a reference conformation. The LP-RMSD is the
    minimum root-mean squared deviation between two sets of points, minimizing
    over both the rotational/translational degrees of freedom AND the label
    correspondences between points in the target and reference conformations.
    This means that it can be used meaningfully with atoms with exchange
    symmetry like, like multiple water molecules.

    Notes
    -----
    The optimization procedure used does necessary not find the rotation/mapping
    for the global minimum of the means squared deviation. It does usually find
    a *pretty good* solution though, using a 3 step procedure. For each frame in
    the target, (1) align the target reference frame using only the
    distinguishable atoms (2) using this fixed rotation matrix, find the optimal
    mapping for the labels in the permute groups (3) using this mapping,
    recompute a new optimal rotation matrix to align the frame, and get the
    final RMSD.

    Parameters
    ----------
    target : md.Trajectory
        For each conformation in this trajectory, compute the RMSD to
        a particular 'reference' conformation in another trajectory
        object.
    reference : md.Trajectory
        The object containing the reference conformation to measure distances
        to.
    frame : int
        The index of the conformation in `reference` to measure
        distances to.
    atom_indices : array_like, or None
        The indices of the atoms to use in the RMSD calculation. If not
        supplied, all atoms will be used.
    permute_groups : list of array_like, or None
        A list of groups of permutable atoms. Each element in permute_groups
        is an array of indices containing atoms whose labels can be mutually
        exchanged within the group. All indices in permute_groups must be a
        subset of atom_indices. The default, None, corresponds to a single group
        containing all the atoms in atom_indices.
    parallel : bool
        Use OpenMP to calculate each of the RMSDs in parallel over multiple
        cores.
    superpose : bool
        Modify the trajectory ``target`` in-place to superpose each frame on
        reference.

    Returns
    -------
    lprmsds : np.ndarray, shape=(target.n_frames,)
        A 1-D numpy array of the linear programming root-mean-square deviations
        from each of the conformations in `target` to the `frame`-th
        conformation in `reference`.
    """
    _validate_shapes(target, reference, frame)
    atom_indices = _validate_atom_indices(atom_indices, target.xyz.shape[1])
    permute_groups = _validate_permute_groups(permute_groups, atom_indices)

    # these are the indices of the permute atoms inside the coords array after
    # selecting only the atom indices.
    # i.e. [xyz[atom_indices][i] for in permute_groups_rel]
    cdef vector[vector[int]] permute_groups_stl
    for pgroup in permute_groups:
         permute_groups_stl.push_back(np.searchsorted(atom_indices, pgroup))
    dis_indices = np.setdiff1d(
            np.arange(len(atom_indices)), np.concatenate(permute_groups_stl))
    # print('dis indices', dis_indices)

    cdef int n_atoms_total = target.xyz.shape[1]
    cdef int superpose_ = superpose
    cdef np.ndarray[ndim=3, dtype=np.float32_t, mode='c'] target_xyz = np.asarray(target.xyz, order='c')
    cdef np.ndarray[ndim=1, dtype=int, mode='c'] atom_indices_ = np.asarray(atom_indices, dtype=np.int32)
    cdef np.ndarray[ndim=1, dtype=int, mode='c'] dis_indices_ =  np.asarray(dis_indices, dtype=np.int32)

    # get the all the coords in atom_indices in target and ref.
    # center them and compute the g values
    cdef int i
    cdef double msd
    cdef float ref_g, target_g
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] target_xyz_frame
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] target_xyz_frame2
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] ref_xyz_frame
    cdef int target_n_frames = target.xyz.shape[0]
    cdef int n_atoms = len(atom_indices)
    ref_xyz_frame = np.array(reference.xyz[frame, atom_indices, :], dtype=np.float32, copy=True)
    target_xyz_frame = np.finfo(np.float32).max * np.ones_like(ref_xyz_frame)
    target_xyz_frame2 = np.empty_like(target_xyz_frame)
    inplace_center_and_trace_atom_major(&ref_xyz_frame[0, 0], &ref_g, 1, n_atoms)


    # get only the distinguishable atoms in target and ref, center them
    # and compute the g values
    cdef float ref_g_dis, target_g_dis
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] target_xyz_frame_dis
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] ref_xyz_frame_dis
    cdef int n_atoms_dis = len(dis_indices)
    target_xyz_frame_dis = np.empty((n_atoms_dis, 3), dtype=np.float32)
    if n_atoms_dis > 0:
        # note: the indexing here is subtpe -- we're using atom_indices[dis_indices] since
        # dis_indices have different semantics from atom_indices
        ref_xyz_frame_dis = np.array(reference.xyz[frame, atom_indices[dis_indices], :], copy=True, order='c')
        inplace_center_and_trace_atom_major(&ref_xyz_frame_dis[0, 0], &ref_g_dis, 1, n_atoms_dis)


    cdef vector[int] mapping
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] rot1 = np.eye(3, dtype=np.float32)
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] rot2 = np.eye(3, dtype=np.float32)
    cdef np.ndarray[ndim=2, dtype=np.float32_t, mode='c'] rot3 = np.eye(3, dtype=np.float32)
    cdef np.ndarray[ndim=1, dtype=np.float32_t, mode='c'] distances = np.zeros(target_n_frames, dtype=np.float32)

    #for i in prange(target_n_frames, nogil=True):
    for i in range(target_n_frames):
        fancy_index2d(&target_xyz[i,0,0], n_atoms_total, 3, <int*> &atom_indices_[0],
                      n_atoms, NULL, 0, &target_xyz_frame[0, 0])
        inplace_center_and_trace_atom_major(&target_xyz_frame[0, 0], &target_g, 1, n_atoms)

        # compute rotation matrix on distinguishable atoms
        if (n_atoms_dis > 0):
            # subsample the target_xyz_frame to get just the distinguishable atoms
            # target_xyz_frame has already been "sampled down" by atom_indices, so we
            # can now apply dis_indices
            fancy_index2d(&target_xyz_frame[0,0], n_atoms, 3, <int*> &dis_indices_[0],
                          n_atoms_dis, NULL, 0, &target_xyz_frame_dis[0, 0])
            inplace_center_and_trace_atom_major(&target_xyz_frame_dis[0,0], &target_g_dis, 1, n_atoms_dis)

            # get `rot1`, the rotation matrix thats optimal for the distinguishable indices
            msd_atom_major(n_atoms_dis, n_atoms_dis,
                &target_xyz_frame_dis[0, 0], &ref_xyz_frame_dis[0, 0],
                target_g_dis, ref_g_dis, 1, &rot1[0, 0])

            # apply rot1 to all the atom_indices
            rot_atom_major(n_atoms, &target_xyz_frame[0, 0], &rot1[0, 0])

        # compute the optimal remapping of the indices
        mapping = euclidean_permutation(&ref_xyz_frame[0, 0], &target_xyz_frame[0, 0], n_atoms, 3, permute_groups_stl)
        # and remap the indices -- i.e.   target_xyz_frame2 = target_xyz_frame[mapping]
        fancy_index2d(&target_xyz_frame[0, 0], n_atoms, 3, &mapping[0], n_atoms, NULL, 0, &target_xyz_frame2[0, 0])

        # msd = ((target_xyz_frame2 - ref_xyz_frame)**2).sum(1).mean(0)
        # then using these remapped indices, compute the rmsd, with a new rotation
        if superpose_:
            msd = msd_atom_major(n_atoms, n_atoms, &target_xyz_frame2[0, 0], &ref_xyz_frame[0, 0], target_g, ref_g, 1, &rot2[0, 0])
            inplace_center_and_trace_atom_major(&target_xyz[i, 0, 0], NULL, 1, n_atoms_total)
            sgemm33(&rot1[0,0], &rot2[0,0], &rot3[0,0])
            rot_atom_major(n_atoms_total, &target_xyz[i, 0, 0], &rot3[0, 0])
        else:
            msd = msd_atom_major(n_atoms, n_atoms, &ref_xyz_frame[0, 0], &target_xyz_frame2[0, 0], target_g, ref_g, 0, NULL)



        distances[i] = sqrtf(msd)

    if superpose_:
        target.xyz = target_xyz

    return distances


def _validate_atom_indices(atom_indices, n_atoms):
    if atom_indices is None:
        atom_indices = np.arange(n_atoms, dtype=int)
    else:
        atom_indices = ensure_type(np.unique(atom_indices), dtype=int,
                                   ndim=1, name='atom_indices', warn_on_cast=False)
        if not np.all((atom_indices >= 0) * (atom_indices < n_atoms) * (atom_indices < n_atoms)):
            raise ValueError("atom_indices must be valid positive indices")
    return atom_indices


def _validate_shapes(target, reference, frame):
    assert (target.xyz.ndim == 3) and (reference.xyz.ndim == 3) and \
           (target.xyz.shape[2]) == 3 and (reference.xyz.shape[2] == 3)
    if not (target.xyz.shape[1]  == reference.xyz.shape[1]):
        raise ValueError("Input trajectories must have same number of atoms. "
                         "found %d and %d." % (target.xyz.shape[1], reference.xyz.shape[1]))
    if frame >= reference.xyz.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: reference has "
                         "only %d frames." % (frame, reference.xyz.shape[0]))


def _validate_permute_groups(permute_groups, atom_indices):
    if permute_groups is None:
        return [atom_indices]

    permute_groups = [ensure_type(np.unique(group), dtype=int, ndim=1,
                         warn_on_cast=False, name='permute_groups[%d]' % i) for \
                         i, group in enumerate(permute_groups)]
    for pgroup in permute_groups:
        if len(np.setdiff1d(pgroup, atom_indices)) > 0:
            raise ValueError('The elements in each permute group must be a subset of atom_indices')

    all_permutable_atoms = np.concatenate(permute_groups)
    if not np.all(np.unique(all_permutable_atoms) == np.sort(all_permutable_atoms)):
        raise ValueError('permute_groups must be mutually disjoint sets')

    return permute_groups


@cython.boundscheck(False)
def _munkres(np.ndarray[np.double_t, ndim=2, mode="c"] A not None):
    """_munkres(A)

    Calculate the minimum cost assignment of a cost matrix, A

    Parameters
    ----------
    A : np.ndarray, dtype=np.double, ndim=2

    Returns
    -------
    assignments : np.ndarray, ndim=2, dtype=int32
        Boolean array with shape equal to the shape of A. assignments[i,j] == 1
        for an assignment, and 0 for a non-assignment

    Examples
    -------
    >>> _munkres(np.array([[7, 4, 3], [6, 8, 5], [9, 4, 4]], dtype=np.double))
    [[0 0 1],
     [1 0 0],
     [0 1 0]]
    """
    cdef int x = A.shape[0]
    cdef int y = A.shape[1]
    cdef np.ndarray[ndim=2, dtype=np.int32_t, mode='c'] rslt

    rslt = np.zeros(shape=(x,y), dtype=np.int32, order='c')
    cdef Munkres* munk = new Munkres()
    munk.solve(<double *> A.data, <int *> rslt.data, x, y)
    del munk

    return rslt
