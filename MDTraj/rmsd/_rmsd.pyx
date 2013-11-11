##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
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


# This code was written by Kyle Beauchamp with some contributes by
# Robert McGibbon

##############################################################################
# Imports
##############################################################################

import cython
import numpy as np

cimport numpy as np
from cpython cimport bool
from cython.parallel cimport prange

np.import_array()

##############################################################################
# External Declarations
##############################################################################

cdef extern float msd_axis_major(int nrealatoms, int npaddedatoms, int rowstride,
                                 float* aT, float* bT, float G_a, float G_b) nogil

cdef extern float msd_atom_major(int nrealatoms, int npaddedatoms,  float* a,
                                 float* b, float G_a, float G_b, int computeRot,
                                 float rot[9]) nogil
cdef extern float rot_msd_atom_major(const int n_real_atoms, const int n_padded_atoms,
                                     const float* a, const float* b, const float rot[9]) nogil

cdef extern from "math.h":
    float sqrtf(float x) nogil


##############################################################################
# Functions
##############################################################################


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleRMSDs_axis_major(
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz1 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz2 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g1 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g2 not None,
int n_atoms,
int frame,
bool parallel=True):
    """getMultipleRMSDs_axis_major(xyz1, xyz2, g1, g2, n_atoms, frame, parallel=True)

    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in axis-major orders

    Parameters
    ----------
    xyz1 : np.ndarray, shape=(n_frames, 3, n_atoms_padded), dtype=float32
        Coordinates of reference frame.
    xyz2 : np.ndarray, shape=(n_frames, 3, n_atoms_padded), dtype=float32
        Coordinates of frames to iterate over
    g1 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz1
    g2 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz2
    n_atoms : int
        The number of atoms in the system.
    frame : int
        Index of the desired reference frame in xyz1.
    parallel : bool, default True
        If True, use OpenMP parallelization.

    Returns
    -------
    rmsds: np.ndarray, shape=(len(xyz2),)
        RMSDS between xyz1[frame] and all of xyz2
    """

    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms_padded = xyz1.shape[2]
    cdef float msd

    assert xyz1.ndim == 3 and xyz2.ndim == 3 and xyz1.shape[1] == 3 and xyz2.shape[1] == 3
    if not ((xyz1.shape[2] % 4 == 0) & (xyz2.shape[2] % 4 == 0)):
        raise ValueError("Input arrays must have third dimension of 4*n, "
                         "found %d and %d." % (xyz1.shape[2], xyz2.shape[2]))
    if frame >= xyz1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz1.shape[0]))

    cdef np.ndarray[dtype=np.float32_t, ndim=1] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd = msd_axis_major(n_atoms, n_atoms_padded, n_atoms_padded,
                                   &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd = msd_axis_major(n_atoms, n_atoms_padded, n_atoms_padded,
                                   &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i])
            distances[i] = sqrtf(msd)

    return distances


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleRMSDs_atom_major(
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz1 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz2 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g1 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g2 not None,
int n_atoms,
int frame,
bool parallel=True):
    """getMultipleRMSDs_atom_major(xyz1, xyz2, g1, g2, n_atoms, frame, parallel=True)

    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in atom-major orders

    Parameters
    ----------
    xyz1 : np.ndarray, shape=(n_frames, n_atoms_padded, 3), dtype=float32
        Coordinates of reference frame.
    xyz2 : np.ndarray, shape=(n_frames, n_atoms_padded, 3), dtype=float32
        Coordinates of frames to iterate over
    g1 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz1
    g2 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz2
    n_atoms : int
        The number of atoms in the system.
    frame : int
        Index of the desired reference frame in xyz1.
    parallel : bool, default True
        If True, use OpenMP parallelization.

    Returns
    -------
    rmsds: np.ndarray, shape=(len(xyz2),)
        RMSDS between xyz1[frame] and all of xyz2
    """

    cdef Py_ssize_t i
    cdef int n_frames = xyz2.shape[0]
    cdef int n_atoms_padded = xyz1.shape[1]
    cdef float msd

    assert xyz1.ndim == 3 and xyz2.ndim == 3 and xyz1.shape[2] == 3 and xyz2.shape[2] == 3
    if not ((xyz1.shape[1] % 4 == 0) & (xyz2.shape[1] % 4 == 0)):
        raise ValueError("Input arrays must have middle dimension of 4*n, "
                         "found %d and %d." % (xyz1.shape[1], xyz2.shape[1]))
    if frame >= xyz1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz1.shape[0]))

    cdef np.ndarray[dtype=np.float32_t, ndim=1] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd = msd_atom_major(n_atoms, n_atoms_padded,  &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i], 0, NULL)
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd = msd_atom_major(n_atoms, n_atoms_padded, &xyz1[frame, 0, 0], &xyz2[i, 0, 0], g1[frame], g2[i], 0, NULL)
            distances[i] = sqrtf(msd)

    return distances


@cython.boundscheck(False)
@cython.wraparound(False)
def getMultipleAlignDisplaceRMSDs_atom_major(
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz_align1 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz_align2 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g_align1 not None,
np.ndarray[np.float32_t, ndim=1, mode="c"] g_align2 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz_displ1 not None,
np.ndarray[np.float32_t, ndim=3, mode="c"] xyz_displ2 not None,
int n_atoms_align,
int n_atoms_displ,
int frame,
bool parallel=True):
    """getMultipleAlignDisplaceRMSDs_atom_major(xyz_align1, xyz_align2, g_align1, g_align2, xyz_displ1, xyz_displ2, n_atoms_align, n_atoms_displ, frame, parallel=True)

    Calculate the RMSD of several frames to a single frame, with the
    coordinates laid out in atom-major orders, where alignments are computed
    by determining the optimal rotation matrix to align `xyz_align1` against
    `xyz_align2`, and then using this rotation matrix displace `xyz_displ1`,
    and computing its mean-squared-displacement with respect to `xyz_displ2`.

    Parameters
    ----------
    xyz_align1 : np.ndarray, shape=(n_frames, n_atoms_align_padded, 3), dtype=float32
        Alignment coordinates of reference frame.
    xyz_align2 : np.ndarray, shape=(n_frames, n_atoms_align_padded, 3), dtype=float32
        Alignment coordinates of frames to iterate over.
    g_align1 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz1
    g_align2 : np.ndarray, shape = (n_frames), dtype=float32
        Pre-calculated G factors (traces) for each frame in xyz2
    xyz_displ1 : np.ndarray, shape=(n_frames, n_atoms_displ_padded, 3), dtype=float32
        Displacement coordinates of reference frame.
    xyz_displ2 : np.ndarray, shape=(n_frames, n_atoms_displ_padded, 3), dtype=float32
        Displacement coordinates of frames to iterate over.
    n_atoms_align : int
        The number of true atoms in the alignment coordinates
    n_atoms_displ : int
        The number of true atoms in the displacement coordinates
    frame : int
        Index of the desired reference frame in xyz_align1 / xyz_displ1.
    parallel : bool, default True
        If True, use OpenMP parallelization.

    Returns
    -------
    rmsds: np.ndarray, shape=(len(xyz2),)
        RMSDs between xyz_displ1[frame] and all of xyz_displ2, after the alignment that maximizes
        the overlap between xyz_align1[frame] and all of xyz_align2
    rotations: np.ndarray, shape=(len(xyz2), 3, 3)
        The rotation matrix for each of the alignments
    """
    assert (xyz_align1.shape[2] == 3) and (xyz_align2.shape[2] == 3)
    assert (xyz_displ1.shape[2] == 3) and (xyz_displ2.shape[2] == 3)
    if not ((xyz_align1.shape[1] % 4 == 0) & (xyz_align2.shape[1] % 4 == 0)):
        raise ValueError("Input arrays must have middle dimension of 4*n, "
                         "found %d and %d." % (xyz_align1.shape[1], xyz_align2.shape[1]))
    if not ((xyz_displ1.shape[1] % 4 == 0) & (xyz_displ2.shape[1] % 4 == 0)):
        raise ValueError("Input arrays must have middle dimension of 4*n, "
                         "found %d and %d." % (xyz_displ1.shape[1], xyz_displ2.shape[1]))
    if not xyz_align1.shape[0] == xyz_displ1.shape[0]:
        raise ValueError("xyz_align1 and xyz_displ1 must contain the same number of frames")
    if not xyz_align2.shape[0] == xyz_displ2.shape[0]:
        raise ValueError("xyz_align2 and xyz_displ2 must contain the same number of frames")
    if frame >= xyz_align1.shape[0]:
        raise ValueError("Cannot calculate RMSD of frame %d: xyz1 has "
                         "only %d frames." % (frame, xyz_align1.shape[0]))

    cdef float msd
    cdef Py_ssize_t i
    cdef int n_frames = xyz_align2.shape[0]
    cdef int n_align_atoms_padded = xyz_align1.shape[1]
    cdef int n_displ_atoms_padded = xyz_displ1.shape[1]

    cdef np.ndarray[ndim=3, dtype=np.float32_t] rot = np.zeros((n_frames, 3, 3), dtype=np.float32)
    cdef np.ndarray[dtype=np.float32_t, ndim=1] distances = np.zeros(n_frames, dtype=np.float32)

    if parallel == True:
        for i in prange(n_frames, nogil=True):
            msd_atom_major(n_atoms_align, n_align_atoms_padded, &xyz_align1[frame, 0, 0], &xyz_align2[i, 0, 0], g_align1[frame], g_align2[i], 1, &rot[i, 0, 0])
            msd = rot_msd_atom_major(n_atoms_displ, n_displ_atoms_padded, &xyz_displ1[frame, 0, 0], &xyz_displ2[i, 0, 0], &rot[i, 0, 0])
            distances[i] = sqrtf(msd)
    else:
        for i in range(n_frames):
            msd_atom_major(n_atoms_align, n_align_atoms_padded, &xyz_align1[frame, 0, 0], &xyz_align2[i, 0, 0], g_align1[frame], g_align2[i], 1, &rot[i, 0, 0])
            msd = rot_msd_atom_major(n_atoms_displ, n_displ_atoms_padded, &xyz_displ1[frame, 0, 0], &xyz_displ2[i, 0, 0], &rot[i, 0, 0])
            distances[i] = sqrtf(msd)

    return distances, rot
