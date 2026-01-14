##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Shehan M. Parmar
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

import numpy as np

__all__ = ["compute_clusters", "compute_cluster_metrics"]


def compute_clusters(traj, residue, criteria, periodic=True):
    """Find molecular clusters based on distance/angle criteria.

    Parameters
    ----------
    traj : md.Trajectory
        Trajectory to analyze.
    residue : str
        Residue name to cluster (e.g., 'HOH').
    criteria : dict or str
        Criteria for connecting molecules.
    periodic : bool, optional, default=True
        Whether to use periodic boundary conditions.

    Returns
    -------
    edges : list of list of tuple
        Per-frame edge lists. Each edge (i, j) uses indices into
        residue_indices (not actual residue indices).
    residue_indices : list of int
        Actual residue indices from topology for the filtered molecules.
        Use residue_indices[i] to get the topology residue index for
        molecule i in the edge list.
    """
    raise NotImplementedError("compute_clusters not yet implemented")


def _build_graph(n_molecules, edges):
    """Build NetworkX graph from molecule count and edge list."""
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_molecules))
    G.add_edges_from(edges)
    return G


def _compute_frame_metrics(graph):
    """Compute cluster sizes and diameters for a single frame."""
    import networkx as nx

    sizes = []
    diameters = []
    for cluster in nx.connected_components(graph):
        sizes.append(len(cluster))
        diameters.append(nx.diameter(graph.subgraph(cluster)))
    return np.array(sizes), np.array(diameters)


def compute_cluster_metrics(edges, n_molecules, return_graphs=False, graph_every=1):
    """Compute cluster sizes and diameters from edge lists.

    Parameters
    ----------
    edges : list of list of tuple
        Edge lists, one list of (i, j) tuples per frame, where i, j 
        are internal indices spanning range(n_molecules). 
    n_molecules : int
        Number of molecules with residues matching the specified criteria.
    return_graphs : bool, optional, default=False
        If True, also return the NetworkX graph objects.
    graph_every : int, optional, default=1
        When return_graphs=True, only return graphs every N frames.

    Returns
    -------
    sizes : list of np.ndarray, length=n_frames
        Cluster sizes for each frame. Each element has shape=(n_clusters,)
        and contains a list of molecule counts per cluster.
    diameters : list of np.ndarray
        Graph diameters per cluster, per frame.
    graphs : list of networkx.Graph, optional
        Only returned if return_graphs=True.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError(
            "NetworkX is required for compute_cluster_metrics. "
            "Install it with: pip install networkx"
        )

    all_sizes = []
    all_diameters = []
    all_graphs = [] if return_graphs else None

    for frame_idx, frame_edges in enumerate(edges):
        G = _build_graph(n_molecules, frame_edges)
        sizes, diameters = _compute_frame_metrics(G)
        all_sizes.append(sizes)
        all_diameters.append(diameters)

        if return_graphs and frame_idx % graph_every == 0:
            all_graphs.append(G)

    if return_graphs:
        return all_sizes, all_diameters, all_graphs
    return all_sizes, all_diameters






