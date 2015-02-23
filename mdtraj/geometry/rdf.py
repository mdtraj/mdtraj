##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
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

import itertools
import re

import numpy as np

from mdtraj.utils import ensure_type
from mdtraj.geometry.distance import compute_distances


__all__ = ['compute_rdf']


def compute_rdf(traj, pair_names=None, r_range=None, bin_width=0.005,
                periodic=True, opt=True):
    """Compute radial distribution functions for pairs in every frame.

    Parameters
    ----------
    traj : Trajectory
    pair_names : array-like, shape=(2,), optional, default=None
        Pair of atom names to consider. Atom names are matched using the same
        regex matching employed by MDTraj's  atom selection DSL. If pair_names
        is not provided, the RDF will be computed between all atoms.
    r_range : array-like, shape=(2,), optional, default=(0.0, 1.0)
        Minimum and maximum radii.
    bin_width : int, optional, default=0.005
        Width of the bins in nanometers.
    opt : bool, default=True
        Use an optimized native library to compute the pair wise distances.

    Returns
    -------
    r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radii values corresponding to the centers of the bins.
    g_r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radial distribution function values at r.
    """
    if not pair_names:
        # all-all
        pairs = np.array(list(itertools.combinations(range(traj.n_atoms), 2)))
    elif len(pair_names) != 2:
        raise ValueError('pair_names must contain two entries if you want to '
                         'calculate the RDF for specific types of atoms.')
    else:
        type_a = [a.index for a in traj.top.atoms
                  if (re.match(pair_names[0], a.name) is not None)]
        type_b = [a.index for a in traj.top.atoms
                  if (re.match(pair_names[1], a.name) is not None)]
        non_matching_atom_names = [pair_names[i] for i, t in enumerate((type_a, type_b))
                                   if len(t) == 0]
        if non_matching_atom_names:
            raise ValueError('Unable to find atoms matching the following '
                             'selection(s): {0}'.format(non_matching_atom_names))
        # TODO: Probably a cleaner way to generate these pairs.
        prod = itertools.product(type_a, type_b)
        pairs = set()
        for a, b in prod:
            if a != b and (b, a) not in pairs:
                pairs.add((a, b))
        pairs = np.array(list(pairs))

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
