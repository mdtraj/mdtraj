##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Ryan DiRisio
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


from collections import defaultdict

import numpy as np

import mdtraj as md
from mdtraj.geometry import compute_distances_core

__all__ = ["pi_stacking"]


def compute_centroid(trajectory, atoms):
    """
    Compute the centroid (center of geometry) of the ring atoms in the trajectory.

    NOTE: Assumes the atoms are not across periodic boundaries.

    Parameters
    ----------
    trajectory : md.Trajectory
        The trajectory to analyze, uses just the xyz coordinates.
    atoms : array_like
        The atom indices involved in the ring.

    Returns
    -------
    centroid : ndarray, shape=[n_frames, 3], dtype=float
        The mean (centroid) of the ring atoms coordinates. MDTraj stores the coordinates
        in nm.
    """
    trj_subset = trajectory.xyz[:, atoms]
    centroid = np.mean(trj_subset, axis=1)
    return centroid


def compute_ring_normal(trajectory, atoms, centroid):
    """
    Calculate the normal vector to the plane defined by the ring atoms and the centroid.

    Parameters
    ----------
    trajectory : md.Trajectory
        The trajectory to analyze, uses just the xyz coordinates.
    atoms : array_like
        The atom indices involved in the ring. Must be at least of length 2.
    centroid : ndarray
        The centroid of the ring atoms.

    Returns
    -------
    norm_vec : np.array, shape=[n_frames, 2, 3], dtype=float

        The normal vector to the plane defined by the ring atoms and the centroid.
        The normal vector is normalized, and the directionality is arbitrarily assigned
        by the order of the atoms list.
    """
    cross = np.cross(
        trajectory.xyz[:, atoms[0]] - centroid,
        trajectory.xyz[:, atoms[1]] - centroid,
    )
    # Normalize vector
    norm_vec = cross / np.linalg.norm(cross, axis=1)[:, np.newaxis]
    return norm_vec


def dot_pdt(v1, v2):
    """
    Compute the dot product of two stacks of vectors.

    ASsumes v1 and v2 are the same shape.

    Parameters
    ----------
    v1 : ndarray
        The first stack of vectors.
    v2 : ndarray
        The second stack of vectors.

    Returns
    -------
    dot : np.array, shape=v1.shape, dtype=float
        The dot product of the two stacks of vectors.
    """
    dot = np.einsum("ij,ij->i", v1, v2)
    return dot


def compute_angles(vec_1, vec_2):
    """
    Compute the angles between the two supplied vectors.

    Assumes the vectors are the same shape.

    Parameters
    ----------
    vec_1 : ndarray
        The first vector to compute the angles between.
    vec_2 : ndarray
        The second vector to compute the angles between.

    Returns
    -------
    angle : np.array, shape=vec_1.shape, dtype=float
        The angles between the two vectors, in radians.
    """
    angle = np.arccos(dot_pdt(vec_1, vec_2) / (np.linalg.norm(vec_1, axis=1) * np.linalg.norm(vec_2, axis=1)))
    return angle


def cap_angle(angle):
    """
    Bring the angle to the range 0-pi/2.

    For this capping, if the angle is beyond pi/2, we are then functionally finding the
    angle where v2 is the negative of the original.

    Parameters
    ----------
    angle : ndarray
        The vector of angles to cap. Angles assumed to be defined between 0 and pi.

    Returns
    -------
    angle :np.array, shape=angle.shape, dtype=float
        The capped angles. Modifies the input array.
    """
    msk = angle > np.pi / 2
    angle[msk] = np.pi - angle[msk]
    return angle


def calculate_intersection_point(
    plane_normal,
    plane_centroid,
    tilted_normal,
    tilted_centroid,
):
    """
    Calculate the point of intersection between the two potentially T-stacked rings.

    Essentially cargo-culted from ProLIF. This implementation may be improved by using
    the actual planes of the rings instead of the normals to start, since then one
    could just solve for the line of intersection, and then project the centroid
    onto that.

    Parameters
    ----------
    plane_normal : ndarray
        The normal vector of one of the rings.
    plane_centroid : ndarray
        The centroid of one of the rings.
    tilted_normal : ndarray
        The normal vector of the other ring.
    tilted_centroid : ndarray
        The centroid of the other ring.

    Returns
    -------
    points : ndarray
        The point of intersection between the two rings. If the two planes are parallel,
        returns array of NaNs for that frame.
    """
    points = np.zeros(plane_normal.shape)
    intersect_direction = np.cross(plane_normal, tilted_normal, axis=1)
    A = np.stack((plane_normal, tilted_normal, intersect_direction), axis=1)
    dets = np.linalg.det(A)
    singular_matrices = np.isclose(dets, 0)
    A[singular_matrices] = np.random.random((3, 3))
    tilted_offset = dot_pdt(tilted_normal, tilted_centroid)
    plane_offset = dot_pdt(plane_normal, plane_centroid)
    d = np.stack(
        (plane_offset, tilted_offset, np.broadcast_to([0.0], plane_offset.shape)),
        axis=1,
    )
    # NOTE: np.linalg.solve changed in v2.0. This is how to do it in 2.0
    intersect_pt = np.linalg.solve(A, d[..., None])[..., 0]
    vec = plane_centroid - intersect_pt
    intersect_direction = intersect_direction / np.linalg.norm(intersect_direction, axis=1)[:, np.newaxis]
    scalar_proj = dot_pdt(intersect_direction, vec)
    non_singular_matrices = ~singular_matrices
    points[non_singular_matrices] = (
        intersect_pt[non_singular_matrices]
        + intersect_direction[non_singular_matrices] * scalar_proj[non_singular_matrices, np.newaxis]
    )
    points[singular_matrices] = np.array([np.nan, np.nan, np.nan])
    return points


def calculate_face_stack_threshold(
    face_plane_angle_range_rad,
    face_normal_to_centroid_angle_range_rad,
    max_face_to_face_centroid_distance,
    plane_angles,
    res_to_lig_angle,
    lig_to_res_angle,
    centroid_dist,
):
    """
    Calculate the threshold for face-to-face pi-stacking interactions.

    Parameters
    ----------
    face_plane_angle_range_rad : tuple of float
        The range of acceptable angles between the two normal vectors defined by the two
        centroids.
    face_normal_to_centroid_angle_range_rad : tuple of float
        The range of acceptable angles between the normal vector and the vector between
        centroids.
        Checks both the ligand and receptor groups.
    max_face_to_face_centroid_distance : float
        The maximum distance between the centroids of the ligand and receptor groups for
        the interaction.
    plane_angles : ndarray
        The angles between the normals of the two rings. In radians.
    res_to_lig_angle : ndarray
        The angles between the normal of the receptor ring and the vector between the
        centroids. In radians.
    lig_to_res_angle : ndarray
        The angles between the normal of the ligand ring and the vector between the
        centroids. In radians.
    centroid_dist : ndarray
        The distance between the centroids of the two rings. In nm.

    Returns
    -------
    face_msk : ndarray
        The mask of frames that meet the threshold for face-to-face pi-stacking.
    """
    face_msk = (
        (face_plane_angle_range_rad[0] <= plane_angles)
        & (plane_angles <= face_plane_angle_range_rad[1])
        & (
            (
                (face_normal_to_centroid_angle_range_rad[0] <= res_to_lig_angle)
                & (res_to_lig_angle <= face_normal_to_centroid_angle_range_rad[1])
            )
            | (
                (face_normal_to_centroid_angle_range_rad[0] <= lig_to_res_angle)
                & (lig_to_res_angle <= face_normal_to_centroid_angle_range_rad[1])
            )
        )
        & (centroid_dist <= max_face_to_face_centroid_distance)
    )
    return face_msk


def calculate_edge_stack_threshold(
    edge_plane_angle_range_rad,
    edge_normal_to_centroid_angle_range_rad,
    max_edge_to_face_centroid_distance,
    edge_intersection_radius,
    plane_angles,
    res_to_lig_angle,
    lig_to_res_angle,
    centroid_dist,
    min_inter_dist,
):
    """
    Calculate the threshold for edge-to-face pi-stacking interactions.

    Parameters
    ----------
    edge_plane_angle_range_rad : tuple of float
        The range of acceptable angles between the two normal vectors defined by the two
        centroids.
    edge_normal_to_centroid_angle_range_rad : tuple of float
        The range of acceptable angles between the normal vector and the vector between
        centroids. Checks both the ligand and receptor groups.
    max_edge_to_face_centroid_distance : float
        The maximum distance between the centroids of the ligand and receptor groups for
        the interaction.
    edge_intersection_radius : float
        The maximum distance between the point of intersection between both rings and
        the opposite ring's centroid.
    plane_angles : ndarray
        The angles between the normals of the two rings. In radians.
    res_to_lig_angle : ndarray
        The angles between the normal of the receptor ring and the vector between the
        centroids. In radians.
    lig_to_res_angle : ndarray
        The angles between the normal of the ligand ring and the vector between the
        centroids. In radians.
    centroid_dist : ndarray
        The distance between the centroids of the two rings. In nm.
    min_inter_dist : ndarray
        The distance between the intersection point and the centroids of the two rings.
        In nm.

    Returns
    -------
    tstack_msk : np.array, shape=[n_frames], dtype=float
        The mask of frames that meet the threshold for edge-to-face pi-stacking.
    """
    tstack_msk = (
        (edge_plane_angle_range_rad[0] <= plane_angles)
        & (plane_angles <= edge_plane_angle_range_rad[1])
        & (
            (
                (edge_normal_to_centroid_angle_range_rad[0] <= res_to_lig_angle)
                & (res_to_lig_angle <= edge_normal_to_centroid_angle_range_rad[1])
            )
            | (
                (edge_normal_to_centroid_angle_range_rad[0] <= lig_to_res_angle)
                & (lig_to_res_angle <= edge_normal_to_centroid_angle_range_rad[1])
            )
        )
        & (centroid_dist <= max_edge_to_face_centroid_distance)
        & (min_inter_dist <= edge_intersection_radius)
    )
    return tstack_msk


def pi_stacking(
    trajectory,
    ligand_aromatic_groups,
    receptor_aromatic_groups,
    ligand_neighbor_cutoff=None,
    max_face_to_face_centroid_distance=5.5,
    face_plane_angle_range=(0.0, 35.0),
    face_normal_to_centroid_angle_range=(0.0, 33.0),
    max_edge_to_face_centroid_distance=0.65,
    edge_plane_angle_range=(50.0, 90.0),
    edge_normal_to_centroid_angle_range=(0.0, 30.0),
    edge_intersection_radius=0.15,
):
    """
    Calculate the pi-stacking interactions based on supplied atom groups.

    This function calculates the pi-stacking interactions between the ligand and
    receptor groups based on the supplied atom indices. The function uses the
    centroid of the groups to determine the distance between them, and the normal
    vector of the groups to determine the angle between them.

    Two types of pi-stacking interactions are considered:

    1. Face-to-face interactions: The two groups are parallel and close to each other.
    2. Edge-to-face interactions: The two groups are not parallel, but the edge of one
       group is close to the face of the other group.

    Both are returned together in the same list.


    Parameters
    ----------
    trajectory : md.Trajectory
        The trajectory to analyze. The ligand and receptor groups map onto the
          trajectory's topology indices.
    ligand_aromatic_groups : list of int
        The atom indices of the groups to be considered aromatic for the ligand.
    receptor_aromatic_groups : list of int
        The atom indices of the groups to be considered aromatic for the receptor.
    ligand_neighbor_cutoff : float, default=None
        The distance cutoff for considering a receptor group for pi-stacking with a
        ligand group. If None, then all pairwise aromatic groups are considered.
        IMPORTANT: When provided, this filtering is only applied once, using the first
        frame of the trajectory, which may miss interactions that form later due to
        molecular movement. It is recommended that you use this parameter only if:
        (1) your system has minimal conformational changes throughout the trajectory, or
        (2) you're only interested in interactions present in the first frame.
    max_face_to_face_centroid_distance : float
        The maximum distance between the centroids of the ligand and receptor groups for
        the interaction.
    face_plane_angle_range : tuple of float
        The range of acceptable angles between the two normal vectors defined by the two
        centroids. In degrees.
    face_normal_to_centroid_angle_range : tuple of float
        The range of acceptable angles between the normal vector and the vector between
        centroids. Checks both the ligand and receptor groups. In degrees.
    max_edge_to_face_centroid_distance : float
        The maximum distance between the centroids of the ligand and receptor groups for
        the interaction.
    edge_plane_angle_range : tuple of float
        The range of acceptable angles between the two normal vectors defined by the two
        centroids. In degrees.
    edge_normal_to_centroid_angle_range : tuple of float
        The range of acceptable angles between the normal vector and the vector between
        centroids. Checks both the ligand and receptor groups. In degrees.
    edge_intersection_radius : float
        The maximum distance between the point of intersection between both rings and
        the opposite ring's centroid.

    Returns
    -------
    stacking_interactions: list, len=n_frames
        A list of lists of tuples, where each tuple is a pair of aromatic groups that
        are stacking in that frame. The order of the tuple goes
        (ligand_group, protein_group). Includes both face-to-face and edge-to-face
        interactions.
    """
    face_plane_angle_range_rad = tuple(np.deg2rad(face_plane_angle_range))
    face_normal_to_centroid_angle_range_rad = tuple(np.deg2rad(face_normal_to_centroid_angle_range))
    edge_plane_angle_range_rad = tuple(np.deg2rad(edge_plane_angle_range))
    edge_normal_to_centroid_angle_range_rad = tuple(np.deg2rad(edge_normal_to_centroid_angle_range))
    # Ensure order of the groups by converting to list.
    ligand_aromatic_groups = list(ligand_aromatic_groups)
    receptor_aromatic_groups = list(receptor_aromatic_groups)
    # ligand_neighbor_groups = {Ligand group idx: receptor group idxs}
    ligand_neighbor_groups = defaultdict(list)
    if ligand_neighbor_cutoff is not None:
        receptor_nbr_atomatic_groups: set[tuple[int, ...]] = set()
        neighbors = md.compute_neighborlist(
            trajectory,
            ligand_neighbor_cutoff,
        )
        for lig_grp in ligand_aromatic_groups:
            # Find all neighbors for each of the atoms in the ligand group, add to set
            lig_neighbors: set[int] = set()
            for atm_idx in lig_grp:
                lig_neighbors.update(neighbors[atm_idx])
            # For each of the receptor groups, check if any of the atoms are a neighbor
            for rec_grp in receptor_aromatic_groups:
                if any(atm in lig_neighbors for atm in rec_grp):
                    ligand_neighbor_groups[lig_grp].append(rec_grp)
                    receptor_nbr_atomatic_groups.add(rec_grp)
    else:
        receptor_nbr_atomatic_groups = set(receptor_aromatic_groups)
        for rec_grp in receptor_aromatic_groups:
            for lig_grp in ligand_aromatic_groups:
                ligand_neighbor_groups[lig_grp].append(rec_grp)
    stacking_interactions = [[] for _ in range(len(trajectory))]
    # If no receptor aromatic neighbors within any of the ligand groups, return
    if len(ligand_neighbor_groups) == 0:
        return stacking_interactions
    # Calculate the centroid and ring normals of each group
    receptor_grp_centroids = {}
    receptor_grp_normals = {}
    for rec_grp in receptor_nbr_atomatic_groups:
        rec_centroid = compute_centroid(trajectory, rec_grp)
        receptor_grp_centroids[rec_grp] = rec_centroid
        receptor_grp_normals[rec_grp] = compute_ring_normal(trajectory, rec_grp, rec_centroid)
    # For each ligand -- receptor pairs, calculate the centroid distance
    for lig_grp, receptor_grps in ligand_neighbor_groups.items():
        lig_centroid = compute_centroid(
            trajectory,
            lig_grp,
        )
        lig_normal = compute_ring_normal(
            trajectory,
            lig_grp,
            lig_centroid,
        )
        # Calculate the plane angles
        for grp in receptor_grps:
            receptor_centroid = receptor_grp_centroids[grp]
            centroid_dist = compute_distances_core(
                np.hstack(
                    (lig_centroid[:, np.newaxis], receptor_centroid[:, np.newaxis]),
                ),
                [[0, 1]],
                trajectory.unitcell_vectors,
            ).squeeze()
            if np.all(
                centroid_dist
                > max(
                    max_face_to_face_centroid_distance,
                    max_edge_to_face_centroid_distance,
                ),
            ):
                # centroid distance too far for all frames, skip
                continue
            # Calculate the angle between the normals
            # Going to assume PBC doesn't apply here since we're just using unit vectors
            res_normal = receptor_grp_normals[grp]
            plane_angles = cap_angle(compute_angles(lig_normal, res_normal))
            # Calcualte the normal vec - centroid vec angles
            res_to_lig_centroid = lig_centroid - receptor_centroid
            res_to_lig_angle = cap_angle(compute_angles(res_normal, res_to_lig_centroid))
            lig_to_res_centroid = receptor_centroid - lig_centroid
            lig_to_res_angle = cap_angle(compute_angles(lig_normal, lig_to_res_centroid))
            lig_res_frames_msk = calculate_face_stack_threshold(
                face_plane_angle_range_rad,
                face_normal_to_centroid_angle_range_rad,
                max_face_to_face_centroid_distance,
                plane_angles,
                res_to_lig_angle,
                lig_to_res_angle,
                centroid_dist,
            )
            for frame in np.where(lig_res_frames_msk)[0]:
                stacking_interactions[frame].append((lig_grp, grp))

            # T-stacks
            intersections = calculate_intersection_point(
                lig_normal,
                lig_centroid,
                res_normal,
                receptor_centroid,
            )
            lig_inter_dist = np.linalg.norm(lig_centroid - intersections, axis=1)
            rec_inter_dist = np.linalg.norm(receptor_centroid - intersections, axis=1)
            min_inter_dist = np.fmin(
                lig_inter_dist,
                rec_inter_dist,
            )
            tstack_msk = calculate_edge_stack_threshold(
                edge_plane_angle_range_rad,
                edge_normal_to_centroid_angle_range_rad,
                max_edge_to_face_centroid_distance,
                edge_intersection_radius,
                plane_angles,
                res_to_lig_angle,
                lig_to_res_angle,
                centroid_dist,
                min_inter_dist,
            )
            for frame in np.where(tstack_msk)[0]:
                stacking_interactions[frame].append((lig_grp, grp))
    return stacking_interactions
