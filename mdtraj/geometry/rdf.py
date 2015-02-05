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

from itertools import combinations

import numpy as np

from mdtraj.utils import ensure_type
from mdtraj.geometry.distance import compute_distances


__all__ = ['compute_rdf']


def distance_pbc(x0, x1, dimensions):
    """Vectorized distance calculationidering minimum image. """
    d = np.abs(x0 - x1)
    d = np.where(d > 0.5 * dimensions, dimensions - d, d)
    return np.sqrt((d ** 2).sum(axis=-1))


def compute_rdf(traj, pair_names=None, periodic=True, r_range=None, n_bins=100):
    """Compute radial distribution functions for pairs in every frame.

    Parameters
    ----------
    traj : Trajectory
    pair_names : array-like, optional
        Pair of atomtypes to consider.
    r_range : array-like, optional
        Minimum and maximum radii.
    n_bins : int, optional
        Number of bins.

    Returns
    -------
    r : ndarray
        Radii values corresponding to bins.
    g_r : ndarray
        Distribution function values at r.
    """
    n_bins = np.uint32(n_bins)
    if not r_range:
        r_range = np.array([0.0, 8.0])
    r_range = np.array(r_range)
    g_r, edges = np.histogram([0], bins=n_bins, range=r_range)
    g_r[0] = 0

    rho = 0
    all_box_lengths = traj.unitcell_lengths
    for frame, xyz in enumerate(traj.xyz):
        print(frame)
        box_lengths = all_box_lengths[frame]
        print(box_lengths, type(box_lengths))
        for i, xyz_i in enumerate(xyz):
            xyz_j = np.vstack([xyz[:i], xyz[i+1:]])
            d = distance_pbc(xyz_i, xyz_j, box_lengths)
            temp_g_r, _ = np.histogram(d, n_bins, r_range)
            g_r += temp_g_r
        rho += (i + 1) / np.sum(box_lengths)

    # normalization
    #n_atoms = traj.n_atoms * traj.n_frames
    #rho = np.sum(n_atoms / traj.unitcell_volumes)
    r = 0.5 * (edges[1:] + edges[:-1])
    V = 4./3. * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    norm = rho * i
    g_r = g_r.astype(np.float32)  # from uint32
    g_r /= norm * V
    return r, g_r


if __name__ == "__main__":
    import pdb

    import mdtraj as md
    from mdtraj.testing import get_fn

    traj = md.load(get_fn('water_216.lammpstrj'), top=get_fn('spc216.pdb'))
    print('loaded')

    r, g_r = compute_rdf(traj)

    import matplotlib.pyplot as plt
    plt.plot(r, g_r)
    plt.show()



