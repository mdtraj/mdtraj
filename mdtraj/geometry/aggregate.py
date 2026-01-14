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

try:
    import networkx as nx
    HAVE_NETWORKX = True
except ImportError:
    HAVE_NETWORKX = False

__all__ = ["compute_aggregates", "compute_aggregate_metrics"]

def _to_name_set(names):
    """Convert atom name(s) to a set, validating input type."""
    if isinstance(names, str):
        return {names}
    if isinstance(names, list):
        return set(names)
    raise TypeError(f"Atom names must be str or list, got {type(names).__name__}")


def _validate_atom_names(expected, atom_indices, topology, source):
    """Raise ValueError if any expected atom names are missing from atom_indices."""
    found = {topology.atom(a).name for a in atom_indices}
    if missing := expected - found:
        raise ValueError(f"Atom names not found in topology: {missing} (from {source})")

def compute_aggregates(
        traj, 
        residue, 
        criteria, 
        periodic=True
):
    """Find molecular aggregates based on distance/angle criteria.

    Parameters
    ----------
    traj : Trajectory
        Trajectory to conduct aggregation analysis.
    residue : str
        Residue name to filter for all aggregates (e.g., 'HOH').
    criteria : dict or list of dict
        Edge criteria for connecting molecules. Can be a single dict or a list
        of dicts for multiple criteria (ALL criteria must be satisfied for an edge).
        Each dict has keys:
        - 'atom_pair' : tuple of (str or list, str or list)
            Atom names for distance check (e.g., ('H', 'O') or (['H1', 'H2'], 'O')).
        - 'distance' : float
            Distance cutoff in nm.
        - 'atom_triplet' : tuple of (str or list, str or list, str or list), optional
            Atom names for angle check as (atom1, vertex, atom2). The angle
            is measured at the vertex atom (e.g., ('O', ['H1', 'H2'], 'O') measures
            the O-H-O angle with H as the vertex). Must have associated `angle.`
        - 'angle' : float, optional
            Minimum angle threshold in degrees. Angles greater than given value
            pass (e.g., 120 means angles > 120° pass). Must have associated 
            `atom_triplet.`

        Example with multiple criteria::

            criteria = [
                {'atom_pair': (['H1', 'H2'], 'O'), 'distance': 0.19},
                {'atom_pair': ('O', 'O'), 'distance': 0.28},
            ]
    periodic : bool, optional, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, distances will be computed based on minimum 
        image convention.

    Returns
    -------
    edges : list of list of tuple
        Edge lists, one list of (i, j) tuples per frame, where i, j
        are internal indices spanning range(n_molecules), where n_molecules
        are the number of (filtered) residues.
    residue_indices : list of int
        Actual residue indices from topology for the filtered molecules.
        residue_indices[i] gets the topology residue index for molecule
        i in the edge list.
    """
    residue_indices = _get_residue_indices(traj.topology, residue)

    if not residue_indices:
        raise ValueError(f"No residues found matching '{residue}'")

    edges = []
    for frame_idx in range(traj.n_frames):
        frame_edges = _find_molecular_edges_multi(traj, frame_idx, residue_indices, criteria, periodic)
        edges.append(frame_edges)
    return edges, residue_indices

def _get_residue_indices(topology, residue_name):
    """Get indices of residues matching the given name."""
    return [r.index for r in topology.residues if r.name == residue_name]

def _find_molecular_edges(
    traj,
    frame_idx,
    criterion,
    periodic,
    topology,
    residue_indices,
    res_idx_to_mol_idx,
    atom_to_res,
):
    """Find all edges (connected molecule pairs) for a single criterion."""
    from mdtraj import compute_angles, compute_neighborlist

    distance_cutoff = criterion['distance']
    atom1_names, atom2_names = map(_to_name_set, criterion['atom_pair'])

    # create atom type sets for convenient criteria lookup
    atoms_type1 = set()
    atoms_type2 = set()
    for res_idx in residue_indices:
        res = topology.residue(res_idx)
        for a in res.atoms:
            if a.name in atom1_names:
                atoms_type1.add(a.index)
            if a.name in atom2_names:
                atoms_type2.add(a.index)

    _validate_atom_names(atom1_names, atoms_type1, topology, "atom_pair[0]")
    _validate_atom_names(atom2_names, atoms_type2, topology, "atom_pair[1]")

    neighbors = compute_neighborlist(traj, cutoff=distance_cutoff, frame=frame_idx, periodic=periodic)

    edges = set()
    # criteria with only distances specified 
    if 'atom_triplet' not in criterion or criterion.get('angle') is None:
        for atom_i in atoms_type1:
            res_i = atom_to_res[atom_i]
            mol_i = res_idx_to_mol_idx[res_i]
            for atom_j in neighbors[atom_i]:
                if atom_j in atoms_type2:
                    res_j = atom_to_res[atom_j]
                    if res_i != res_j:
                        mol_j = res_idx_to_mol_idx[res_j]
                        edge = (mol_i, mol_j) if mol_i < mol_j else (mol_j, mol_i)
                        edges.add(edge)
        return edges

    triplet_names = criterion['atom_triplet']
    angle_cutoff = criterion['angle']
    triplet_first, triplet_vertex, triplet_third = map(_to_name_set, triplet_names)

    # create atom type sets for each atom in triplet
    triplet_first_atoms = {}
    atoms_triplet_third = set()
    for res_idx in residue_indices:
        res = topology.residue(res_idx)
        triplet_first_atoms[res_idx] = []
        for a in res.atoms:
            if a.name in triplet_first:
                triplet_first_atoms[res_idx].append(a.index)
            if a.name in triplet_third:
                atoms_triplet_third.add(a.index)

    all_triplet_first = [a for atoms in triplet_first_atoms.values() for a in atoms]
    _validate_atom_names(triplet_first, all_triplet_first, topology, "atom_triplet[0]")

    atoms_for_angle = {a for a in atoms_type1 if topology.atom(a).name in triplet_vertex}
    _validate_atom_names(triplet_vertex, atoms_for_angle, topology, "atom_triplet[1]")
    _validate_atom_names(triplet_third, atoms_triplet_third, topology, "atom_triplet[2]")

    # Collect atom pairs that pass distance cutoff
    passing_atom_pairs = []
    for atom_i in atoms_for_angle:
        res_i = atom_to_res[atom_i]
        mol_i = res_idx_to_mol_idx[res_i]
        for atom_j in neighbors[atom_i]:
            if atom_j in atoms_triplet_third:
                res_j = atom_to_res[atom_j]
                if res_i != res_j:
                    mol_j = res_idx_to_mol_idx[res_j]
                    passing_atom_pairs.append((atom_i, atom_j, res_i, mol_i, mol_j))

    if not passing_atom_pairs:
        return edges

    # create list of atom triplets and molecular pairs (that are within distance cutoff)  
    angle_triplets = []
    triplet_mol_pairs = []
    for atom_vertex, atom_acceptor, res_vertex, mol_i, mol_j in passing_atom_pairs:
        for atom_donor in triplet_first_atoms[res_vertex]:
            angle_triplets.append([atom_donor, atom_vertex, atom_acceptor])
            edge = (mol_i, mol_j) if mol_i < mol_j else (mol_j, mol_i)
            triplet_mol_pairs.append(edge)

    if not angle_triplets:
        return edges

    # Vectorized angle calculation over pairs within distance cutoff
    angle_triplets_arr = np.array(angle_triplets)
    angles = compute_angles(traj[frame_idx], angle_triplets_arr, periodic=periodic)[0]
    angles_deg = np.degrees(angles)

    # Add edges that pass angle criterion
    angles_passed = angles_deg > angle_cutoff
    for i, passed in enumerate(angles_passed):
        if passed:
            edges.add(triplet_mol_pairs[i])

    return edges


def _find_molecular_edges_multi(traj, frame_idx, residue_indices, criteria, periodic):
    """Find edges that satisfy ALL criteria (intersection)."""
    criteria_list = [criteria] if isinstance(criteria, dict) else criteria
    topology = traj.topology

    # Map residue idx to molecule idx, and atom idx to residue idx
    res_idx_to_mol_idx = {res_idx: i for i, res_idx in enumerate(residue_indices)}
    atom_to_res = {a.index: res_idx for res_idx in residue_indices
        for a in topology.residue(res_idx).atoms}

    # Find intersection of all criteria
    all_edges = None
    for criterion in criteria_list:
        edges = _find_molecular_edges(
            traj, frame_idx, criterion, periodic,
            topology, residue_indices, res_idx_to_mol_idx, atom_to_res,
        )
        if all_edges is None:
            all_edges = edges
        else:
            all_edges = all_edges & edges

    return list(all_edges if all_edges else set())


def _build_graph(n_molecules, edges):
    """Build NetworkX graph from molecule count and edge list."""
    G = nx.Graph()
    G.add_nodes_from(range(n_molecules))
    G.add_edges_from(edges)
    return G


def _compute_frame_metrics(graph, approx=False):
    """Compute aggregate sizes and diameters for a single frame."""
    if approx:
        from networkx.algorithms.approximation import diameter as approx_diameter

    sizes = []
    diameters = []
    for component in nx.connected_components(graph):
        sizes.append(len(component))
        subgraph = graph.subgraph(component)
        if approx:
            diameters.append(approx_diameter(subgraph))
        else:
            diameters.append(nx.diameter(subgraph))
    return np.array(sizes), np.array(diameters)


def compute_aggregate_metrics(edges, n_molecules, approx=False, return_graphs=False, graph_every=1):
    """Compute aggregate sizes and diameters from edge lists.

    Parameters
    ----------
    edges : list of list of tuple
        Edge lists, one list of (i, j) tuples per frame, where i, j
        are internal indices spanning range(n_molecules).
    n_molecules : int
        Number of molecules with residues matching the specified criteria.
    approx : bool, optional, default=False
        If True, use approximate diameter algorithm at ~100x the speed
        (returns lower bound on diameter).
    return_graphs : bool, optional, default=False
        If True, also return the NetworkX graph objects.
    graph_every : int, optional, default=1
        When return_graphs=True, only return graphs every N frames.

    Returns
    -------
    sizes : list of np.ndarray, length=n_frames
        Aggregate sizes for each frame. Each element has shape=(n_aggregates,)
        and contains the molecule count per aggregate.
    diameters : list of np.ndarray
        Graph diameters per aggregate, per frame (lower bound if approx=True).
    graphs : list of networkx.Graph, optional
        Only returned if return_graphs=True.
    """
    if not HAVE_NETWORKX:
        raise ImportError(
            "NetworkX is required for compute_aggregate_metrics. "
            "Install it with: pip install networkx"
        )

    all_sizes = []
    all_diameters = []
    all_graphs = [] if return_graphs else None

    for frame_idx, frame_edges in enumerate(edges):
        G = _build_graph(n_molecules, frame_edges)
        sizes, diameters = _compute_frame_metrics(G, approx=approx)
        all_sizes.append(sizes)
        all_diameters.append(diameters)

        if return_graphs and frame_idx % graph_every == 0:
            all_graphs.append(G)

    if return_graphs:
        return all_sizes, all_diameters, all_graphs
    return all_sizes, all_diameters






