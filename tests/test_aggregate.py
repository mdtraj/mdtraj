
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
from mdtraj.geometry.aggregate import compute_aggregates
from mdtraj.testing import eq

try:
    import networkx as nx
    HAVE_NETWORKX = True
except ImportError:
    HAVE_NETWORKX = False
if HAVE_NETWORKX:
    from mdtraj.geometry.aggregate import compute_aggregate_metrics

needs_networkx = pytest.mark.skipif(not HAVE_NETWORKX, reason="needs networkx")


def _hbond_criteria(distance, angle=None):
    """Build H-bond criteria dict for water molecules."""
    criteria = {'atom_pair': (['H1', 'H2'], 'O'), 'distance': distance}
    if angle is not None:
        criteria['atom_triplet'] = ('O', ['H1', 'H2'], 'O')
        criteria['angle'] = angle
    return criteria


def _get_aggregates(traj, criteria, periodic=False):
    """Compute aggregates and return (n_aggregates, sizes, n_edges)."""
    edges, residue_indices = compute_aggregates(traj, residue='HOH', criteria=criteria, periodic=periodic)
    sizes, _ = compute_aggregate_metrics(edges, len(residue_indices))
    return len(sizes[0]), tuple(sizes[0]), len(edges[0])


@needs_networkx
def test_aggregate_metrics():
    """Test aggregate sizes and diameters on synthetic graph data."""
    n_molecules = 5
    edges = [[(0, 1), (0, 2), (0, 3), (2, 3), (2, 4), (3, 4)]]

    sizes, diameters = compute_aggregate_metrics(edges, n_molecules)

    eq(len(sizes), 1)
    eq(sizes[0], np.array([5]))
    eq(len(diameters), 1)
    eq(diameters[0], np.array([3]))


@needs_networkx
def test_compute_aggregates_ice(get_fn):
    """Test compute_aggregates with ice structure."""
    t = md.load(get_fn("ice_1c_1x1x1.pdb"))

    # Distance-only
    eq((1, (8,), 7), _get_aggregates(t, _hbond_criteria(0.19)))
    eq((8, (1, 1, 1, 1, 1, 1, 1, 1), 0), _get_aggregates(t, _hbond_criteria(0.17)))

    # With angle criterion
    eq((1, (8,), 7), _get_aggregates(t, _hbond_criteria(0.19, 150.0)))
    eq((8, (1, 1, 1, 1, 1, 1, 1, 1), 0), _get_aggregates(t, _hbond_criteria(0.19, 179.0)))

    # PBC: same aggregate sizes, more edges
    eq((1, (8,), 15), _get_aggregates(t, _hbond_criteria(0.19), periodic=True))


@needs_networkx
def test_compute_aggregates_multi_criteria(get_fn):
    """Test compute_aggregates with multiple criteria (H-O and O-O distances)."""
    t = md.load(get_fn("ice_1c_1x1x1.pdb"))

    def _multi_criteria(ho_distance, oo_distance):
        return [
            {'atom_pair': (['H1', 'H2'], 'O'), 'distance': ho_distance},
            {'atom_pair': ('O', 'O'), 'distance': oo_distance},
        ]

    eq((8, (1, 1, 1, 1, 1, 1, 1, 1), 0), _get_aggregates(t, _multi_criteria(0.17, 0.28)))
    eq((8, (1, 1, 1, 1, 1, 1, 1, 1), 0), _get_aggregates(t, _multi_criteria(0.19, 0.26)))
    eq((1, (8,), 7), _get_aggregates(t, _multi_criteria(0.19, 0.28)))

