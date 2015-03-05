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

from itertools import combinations, chain, product
import re

import numpy as np

from mdtraj.utils import ensure_type
from mdtraj.geometry.distance import compute_distances


__all__ = ['compute_rdf', 'generate_unique_pairs']


def compute_rdf(traj, pairs=None, r_range=None, bin_width=0.005,
                periodic=True, opt=True):
    """Compute radial distribution functions for pairs in every frame.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute radial distribution function in.
    pairs : array-like, shape=(n_pairs, 2), dtype=int, optional, default=None
        Each row gives the indices of two atoms.
    r_range : array-like, shape=(2,), optional, default=(0.0, 1.0)
        Minimum and maximum radii.
    bin_width : int, optional, default=0.005
        Width of the bins in nanometers.
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
    generate_unique pairs

    """
    if not r_range:
        r_range = np.array([0.0, 1.0])
    r_range = ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
                          shape=(2,), warn_on_cast=False)
    bins = np.arange(r_range[0], r_range[1], bin_width)

    distances = compute_distances(traj, pairs, periodic=periodic, opt=opt)
    g_r, edges = np.histogram(distances, bins=bins)
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


def generate_unique_pairs(traj=None, pair_names=None, a_indices=None,
                          b_indices=None):
    """Generate unique pairs of atom indices.

    If only a trajectory is provided, unique pairs between all atoms will be
    returned.

    If only pair_names is provided, unique pairs between all atoms

    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute radial distribution function in.
    pair_names : array-like, shape=(2,), dtype=str, optional, default=None
         Pair of atom names to consider. Atom names are matched using the same
         regex matching employed by MDTraj's atom selection DSL.
    a_indices : array-like, shape=(1, )
    Returns
    -------

    """
    _gave_indices = a_indices is not None and b_indices is not None
    _gave_pair_names = pair_names is not None

    # Provided nothing: generate all-all pairs.
    if not (_gave_pair_names or _gave_indices):
        pairs = np.fromiter(chain.from_iterable(combinations(range(traj.n_atoms), 2)),
                             dtype=np.int32, count=traj.n_atoms * (traj.n_atoms - 1))
        pairs = np.vstack((pairs[::2], pairs[1::2])).T
        return pairs

    assert _gave_pair_names != _gave_indices, (
        'Must provide either pair_names or both a_indices and b_indices.')

    # Find the indices for the provided pair_names.
    if _gave_pair_names:
        pair_names = ensure_type(pair_names, dtype=str, ndim=1, name='pair_names', length=2,
                                 warn_on_cast=False)
        a_indices = [a.index for a in traj.top.atoms
                  if (re.match(pair_names[0], a.name) is not None)]
        b_indices = [a.index for a in traj.top.atoms
                  if (re.match(pair_names[1], a.name) is not None)]
        _non_matching_atom_names = [pair_names[i] for i, t in enumerate((a_indices, b_indices))
                                   if len(t) == 0]
        if _non_matching_atom_names:
            raise ValueError('Unable to find atoms matching the following '
                             'selection(s): {0}'.format(_non_matching_atom_names))

    # Create the pairs from the indices.
    a_indices = ensure_type(a_indices, dtype=np.int32, ndim=1, name='a_indices', warn_on_cast=False)
    b_indices = ensure_type(b_indices, dtype=np.int32, ndim=1, name='b_indices', warn_on_cast=False)
    pairs = np.fromiter(
        chain.from_iterable(
            (a, b) if a > b else (b, a)
            for (a, b) in product(a_indices, b_indices) if a != b),
        dtype=np.int32)
    pairs = np.vstack((pairs[::2], pairs[1::2])).T
    return pairs
