
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
#
# Authors: Shehan M. Parmar
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

import numpy as np
import pytest

import mdtraj as md
from mdtraj.geometry.cluster import compute_cluster_sizes, compute_clusters
from mdtraj.testing import eq

try:
    import networkx as nx
    HAVE_NETWORKX = True
except ImportError:
    HAVE_NETWORKX = False
if HAVE_NETWORKX:
    from mdtraj.geometry.cluster import compute_cluster_diameters

needs_networkx = pytest.mark.skipif(not HAVE_NETWORKX, reason="needs networkx")

def test_cluster_sizes():
    """Test output of `compute_cluster_sizes` via synthetic data."""
    # Two frames: homogeneous (all one cluster) and heterogeneous (3 clusters)
    clusters = [
        np.array([0, 0, 0, 0, 0, 0]),  
        np.array([0, 0, 1, 1, 1, 2]),  
    ]

    sizes = compute_cluster_sizes(clusters)
    eq(len(sizes), 2)
    eq(sizes[0], np.array([6])) # Frame 0: all 6 in cluster 0
    eq(sizes[1], np.array([2, 3, 1])) # Frame 1: sizes 2, 3, 1

@needs_networkx
def test_cluster_diameters():
    """Test diameter calculation on synthetic graph data."""
    cluster_ids = [np.array([0, 0, 0, 0, 0])]  # 5 nodes, all one cluster
    edges = [[(0, 1), (0, 2), (0, 3), (2, 3), (2, 4), (3, 4)]]

    diameters = compute_cluster_diameters(cluster_ids, edges)
    eq(len(diameters), 1)  # one frame
    assert diameters[0][0] == 3

def test_compute_clusters_ice(get_fn):
    """Test compute_clusters with ice structure."""
    t = md.load(get_fn("ice_1c_1x1x1.pdb"))

    def _hbond_criteria(distance, angle=None):
        criteria = {'atom_pair': ('H', 'O'), 'distance': distance}
        if angle is not None:
            criteria['atom_triplet'] = ('O', 'H', 'O')
            criteria['angle'] = angle
        return criteria

    def _get_clusters(criteria):
        cluster_ids, _ = compute_clusters(t, residue='HOH', criteria=criteria, periodic=False)
        sizes = compute_cluster_sizes(cluster_ids)[0]
        return len(sizes), sizes

    eq((1, np.array([8])), _get_clusters(_hbond_criteria(0.19)))
    eq((8, np.ones(8, dtype=int)), _get_clusters(_hbond_criteria(0.17)))

    eq((1, np.array([8])), _get_clusters(_hbond_criteria(0.19, 30.0)))
    eq((8, np.ones(8, dtype=int)), _get_clusters(_hbond_criteria(0.19, 1.0)))

    eq((1, np.array([8])), _get_clusters('baker_hubbard'))
    eq((1, np.array([8])), _get_clusters('wernet_nilsson'))

