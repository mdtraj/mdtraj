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

__all__ = ["compute_cluster_sizes", "compute_cluster_diameters"]

def compute_cluster_sizes(cluster_ids):
    """Compute size (number of molecules) of each unique cluster.

    Parameters
    ----------
    cluster_ids : list of np.ndarray, shape=(n_frames, n_molecules)
        The integer cluster `ID` assigned to each molecule for each frame. 

    Returns
    -------
    sizes : list of np.ndarray, length=n_frames
        Cluster sizes for each frame. Each element has shape=(n_clusters,)
        and contains molecule count per cluster, sourted by cluster ID. 
    """
    sizes = []
    for frame_clusters in cluster_ids:
        _, counts = np.unique(frame_clusters, return_counts=True)
        sizes.append(counts)
    return sizes


def _build_graph(n_nodes, edges):
    """Build NetworkX graph from node count and edge list.

    Parameters
    ----------
    n_nodes : int
        Number of nodes (molecules).
    edges : list of tuple
        List of (i, j) tuples indicating connected nodes.

    Returns
    -------
    G : networkx.Graph
        Graph with nodes 0 to n_nodes-1 and edges added.
    """
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(range(n_nodes))
    G.add_edges_from(edges)
    return G


def _compute_frame_diameters(graph):
    """Compute cluster diameters for a graph.

    Parameters
    ----------
    graph : networkx.Graph
        Graph with nodes and edges.

    Returns
    -------
    diameters : list
        List of graph diameters, one per connected component.
    """
    import networkx as nx

    diameters = []
    for cluster in nx.connected_components(graph):
        diameters.append(nx.diameter(graph.subgraph(cluster)))
    return diameters


def compute_cluster_diameters(cluster_ids, edges, return_graphs=False, graph_every=1):
    """Compute the graph diameter of each cluster.

    The graph diameter is the longest shortest path between any two (molecular) nodes.

    Parameters
    ----------
    cluster_ids : list of np.ndarray, shape=(n_frames, n_molecules)
        The integer cluster `ID` assigned to each molecule for each frame.
    edges : list of list of tuple, shape=(n_frames, n_edges, 2)
        Edge lists, one list of (i, j) tuples per frame.
    return_graphs : bool, optional, default=False
        If True, also return the NetworkX graph objects.
    graph_every : int, optional, default=1
        When return_graphs=True, only return graphs every N frames.

    Returns
    -------
    diameters : list of np.ndarray
        Graph diameters per connected component, per frame.
    graphs : list of networkx.Graph, optional
        Only returned if return_graphs=True.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError(
            "NetworkX is required for compute_cluster_diameters. "
            "Install it with: pip install networkx"
        )

    all_diameters = []
    all_graphs = [] if return_graphs else None

    for frame_idx in range(len(edges)):
        G = _build_graph(len(cluster_ids[frame_idx]), edges[frame_idx])
        diameters = _compute_frame_diameters(G)
        all_diameters.append(np.array(diameters))

        if return_graphs and frame_idx % graph_every == 0:
            all_graphs.append(G)

    if return_graphs:
        return all_diameters, all_graphs
    return all_diameters






