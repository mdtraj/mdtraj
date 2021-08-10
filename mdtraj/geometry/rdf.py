##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2015 Stanford University and the Authors
#
# Authors: Christoph Klein
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

from __future__ import print_function, division

import numpy as np

from mdtraj.utils import ensure_type
from mdtraj.geometry.distance import compute_distances, compute_distances_t

__all__ = ['compute_rdf', 'compute_rdf_t']


def compute_rdf(traj, pairs, r_range=None, bin_width=0.005, n_bins=None,
                periodic=True, opt=True):
    """Compute radial distribution functions for pairs in every frame.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute radial distribution function in.
    pairs : array-like, shape=(n_pairs, 2), dtype=int
        Each row gives the indices of two atoms.
    r_range : array-like, shape=(2,), optional, default=(0.0, 1.0)
        Minimum and maximum radii.
    bin_width : float, optional, default=0.005
        Width of the bins in nanometers.
    n_bins : int, optional, default=None
        The number of bins. If specified, this will override the `bin_width`
         parameter.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.
    opt : bool, default=True
        Use an optimized native library to compute the pair wise distances.

    Returns
    -------
    r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radii values corresponding to the centers of the bins.
    g_r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radial distribution function values at r.

    See also
    --------
    Topology.select_pairs

    """
    if r_range is None:
        r_range = np.array([0.0, 1.0])
    r_range = ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
                          shape=(2,), warn_on_cast=False)
    if n_bins is not None:
        n_bins = int(n_bins)
        if n_bins <= 0:
            raise ValueError('`n_bins` must be a positive integer')
    else:
        n_bins = int((r_range[1] - r_range[0]) / bin_width)

    distances = compute_distances(traj, pairs, periodic=periodic, opt=opt)
    g_r, edges = np.histogram(distances, range=r_range, bins=n_bins)
    r = 0.5 * (edges[1:] + edges[:-1])

    # Normalize by volume of the spherical shell.
    # See discussion https://github.com/mdtraj/mdtraj/pull/724. There might be
    # a less biased way to accomplish this. The conclusion was that this could
    # be interesting to try, but is likely not hugely consequential. This method
    # of doing the calculations matches the implementation in other packages like
    # AmberTools' cpptraj and gromacs g_rdf.
    V = (4 / 3) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    norm = len(pairs) * np.sum(1.0 / traj.unitcell_volumes) * V
    g_r = g_r.astype(np.float64) / norm  # From int64.
    return r, g_r


def compute_rdf_t(traj, pairs, times, period_length=None, r_range=None, 
                  bin_width=0.005, n_bins=None, self_correlation=True, 
                  periodic=True, n_concurrent_pairs = 100000, opt=True):
    """Compute time-dependent radial distribution functions, g(r, t).
    The time-dependent radial distribution function is calculated between pairs of time points.
    For example, g(r, 0) is equal to the time-independent radial distribution function, g(r).
    Please see https://doi.org/10.1103/PhysRev.95.249 for further reference.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute time-dependent radial distribution function on.
    pairs : array-like, shape=(n_pairs, 2), dtype=int
        Each row gives the indices of two atoms.
    times : array-like, shape=(any, 2), dtype=int
        Each row gives the indices of two frames.
    period_length : int, optional, default=None
        The length of each chunk of frames to consider when time-averaging
    r_range : array-like, shape=(2,), optional, default=(0.0, 1.0)
        Minimum and maximum radii.
    bin_width : float, optional, default=0.005
        Width of the bins in nanometers.
    n_bins : int, optional, default=None
        The number of bins. If specified, this will override the `bin_width`
         parameter.
    self_correlation : bool, default=True
        Whether or not to include the self-correlation, the case of i=j
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.
    n_concurrent_pairs : int, default=100000
        Number of atom pairs analyzed at a time.
    opt : bool, default=True
        Use an optimized native library to compute the pair wise distances.

    Returns
    -------
    r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radii values corresponding to the centers of the bins.
    g_r_t : np.ndarray, shape=(len(times), np.diff(r_range) / bin_width - 1), dtype=float
        Radial distribution function values at r for each time pair.

    See also
    --------
    Topology.select_pairs

    """
    if r_range is None:
        r_range = np.array([0.0, 1.0])
    r_range = ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
                          shape=(2,), warn_on_cast=False)
    if n_bins is not None:
        n_bins = int(n_bins)
        if n_bins <= 0:
            raise ValueError('`n_bins` must be a positive integer')
    else:
        n_bins = int((r_range[1] - r_range[0]) / bin_width)

    if period_length is None:
        period_length = traj.n_frames

    # Add self pairs to `pairs`
    if self_correlation:
        pairs_set = np.unique(pairs)
        pairs = np.vstack([np.vstack([pairs_set, pairs_set]).T, pairs])

    n_small_chunks = np.ceil(len(pairs)/n_concurrent_pairs).astype("int")
    g_r_t = np.zeros((n_small_chunks, len(times), n_bins))
    weights = np.zeros(n_small_chunks)

    # Splits pairs into smaller chunks so that frame_distances is not excessively large
    for i in range(n_small_chunks):
        temp_g_r_t = np.zeros((len(times), n_bins))
        pairs_set = pairs[i*n_concurrent_pairs:(i+1)*n_concurrent_pairs]
        weights[i] = len(pairs_set)/n_concurrent_pairs
        
        # Returns shape (len(times), len(pairs_set))
        frame_distances = compute_distances_t(traj, pairs_set, times, periodic=periodic, opt=opt)

        for n, distances in enumerate(frame_distances):
            tmp, edges = np.histogram(distances, range=r_range, bins=n_bins)
            temp_g_r_t[n, :] += tmp
        r = 0.5 * (edges[1:] + edges[:-1])

        # Normalize by volume of the spherical shell (see above)
        V = (4 / 3) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
        norm = len(pairs_set) / (period_length) * np.sum(1.0 / traj.unitcell_volumes) * V

        temp_g_r_t = temp_g_r_t.astype(np.float64) / norm  # From int64.
        g_r_t[i] = temp_g_r_t

    g_r_t_final = np.average(g_r_t, axis=0, weights=weights)
    
    return r, g_r_t_final
