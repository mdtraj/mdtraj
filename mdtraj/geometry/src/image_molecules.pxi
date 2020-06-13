
# This file is lovingly included from _geometry.pyx

import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libc.stdint cimport int32_t, int64_t

cdef extern from "math_patch.h" nogil:
    float roundf(float x)
    float floorf(float x)

cdef void make_whole(float[:,::1] frame_positions,
                float[:,::1] frame_unitcell_vectors,
                int32_t[:,:] sorted_bonds) nogil:
    # Fix each molecule to ensure the periodic boundary conditions are not
    # splitting it into pieces.
    cdef int atom1, atom2, j, k
    cdef float offset[3]
    cdef float delta[3]
    for j in range(sorted_bonds.shape[0]):
        atom1 = sorted_bonds[j, 0]
        atom2 = sorted_bonds[j, 1]
        for k in range(3):
            delta[k] = frame_positions[atom2, k]  - frame_positions[atom1, k]
            offset[k] = frame_unitcell_vectors[2, k]*roundf(delta[2]/frame_unitcell_vectors[2,2])
            offset[k] += frame_unitcell_vectors[1, k]*roundf((delta[1]-offset[1])/frame_unitcell_vectors[1,1])
            offset[k] += frame_unitcell_vectors[0, k]*roundf((delta[0]-offset[0])/frame_unitcell_vectors[0,0])
            frame_positions[atom2, k] = frame_positions[atom2, k] - offset[k]

cdef void anchor_dists(float[:,::1] frame_positions,
                  float[:,::1] frame_unitcell_vectors,
                  vector[int[:]] anchor_molecules,
                  float[:,:] anchor_dist,
                  int[:,:,:] anchor_nearest_atoms,
                  int num_anchors) nogil:
    cdef int ca1, ca2, mol1, mol2
    cdef float cdist
    cdef int[:] atoms1, atoms2
    for mol1 in range(num_anchors):
        atoms1 = anchor_molecules[mol1]
        for mol2 in range(mol1):
            atoms2 = anchor_molecules[mol2]
            find_closest_contact(&frame_positions[0,0], &atoms1[0], &atoms2[0], atoms1.shape[0], atoms2.shape[0], &frame_unitcell_vectors[0,0], &ca1, &ca2, &cdist)
            anchor_dist[mol1, mol2] = cdist
            anchor_dist[mol2, mol1] = cdist
            anchor_nearest_atoms[mol1, mol2, 0] = ca1
            anchor_nearest_atoms[mol1, mol2, 1] = ca2
            anchor_nearest_atoms[mol2, mol1, 0] = ca1
            anchor_nearest_atoms[mol2, mol1, 1] = ca2

cdef void wrap_mols(float[:,::1] frame_positions,
                    float[:,::1] frame_unitcell_vectors,
                    float[:] center,
                    vector[int[:]] other_molecules) nogil:
    # Loop over all molecules, apply the correct offset (so that anchor
    # molecules will end up centered in the periodic box), and then wrap
    # the molecule into the box.
    cdef int i, j, k
    cdef float offset[3]
    for k in range(3):
        offset[k] = 0.5*frame_unitcell_vectors[k,k] - center[k]
    for j in range(frame_positions.shape[0]):
        for k in range(3):
            frame_positions[j, k] += offset[k]

    cdef int[:] mol
    cdef float mol_center[3]
    cdef float mol_offset[3]

    cdef int n
    for i in range(other_molecules.size()):
        mol = other_molecules[i]
        n = mol.shape[0]
        for k in range(3):
            mol_center[k] = 0
        for j in range(n):
            for k in range(3):
                mol_center[k] += frame_positions[mol[j], k]
        for k in range(3):
            mol_center[k] /= n
            mol_offset[k] = mol_center[k]
        for k in range(3):
            mol_offset[k] = mol_center[k] - frame_unitcell_vectors[2, k]*floorf(mol_offset[2]/frame_unitcell_vectors[2,2])
            mol_offset[k] -= frame_unitcell_vectors[1, k]*floorf(mol_offset[1]/frame_unitcell_vectors[1,1])
            mol_offset[k] -= frame_unitcell_vectors[0, k]*floorf(mol_offset[0]/frame_unitcell_vectors[0,0])

        for j in range(n):
            for k in range(3):
                frame_positions[mol[j], k] += mol_offset[k]-mol_center[k]


cdef void image_frame(frame_positions,
                 frame_unitcell_vectors,
                 anchor_molecules,
                 vector[int[:]] other_molecules,
                 sorted_bonds):

    if sorted_bonds is not None:
        make_whole(frame_positions, frame_unitcell_vectors, sorted_bonds)

    # Compute the distance between each pair of anchor molecules in this frame.
    cdef int num_anchors = len(anchor_molecules)
    anchor_dist = np.zeros((num_anchors, num_anchors), dtype=np.float32)
    cdef int[:,:,:] anchor_nearest_atoms = np.zeros((num_anchors, num_anchors, 2), dtype=np.int32)
    anchor_dists(frame_positions, frame_unitcell_vectors, anchor_molecules, anchor_dist, anchor_nearest_atoms, num_anchors)


    # Start by taking the largest molecule as our first anchor.
    used_anchors = np.arange(1)
    available_anchors = np.arange(1, num_anchors)
    min_anchor_dist = anchor_dist[0, :]

    # Add in anchors one at a time, always taking the one that is nearest an existing anchor.
    offset = np.zeros(3, dtype=np.float32)
    cdef int next_anchor, next_index, nearest_to
    while len(available_anchors) > 0:
        next_index = np.argmin(min_anchor_dist[available_anchors])
        next_anchor = available_anchors[next_index]

        # Find which existing anchor it's closest to, and choose the periodic
        # copy that minimizes the distance to that anchor.
        nearest_to = used_anchors[np.argmin(anchor_dist[next_anchor, used_anchors])]
        a1, a2 = anchor_nearest_atoms[next_anchor, nearest_to]
        if a1 in anchor_molecules[next_anchor]:
            a2, a1 = a1, a2
        delta = frame_positions[a2]-frame_positions[a1]
        offset = frame_unitcell_vectors[2]*np.round(delta[2]/frame_unitcell_vectors[2,2])
        offset += frame_unitcell_vectors[1]*np.round((delta[1]-offset[1])/frame_unitcell_vectors[1,1])
        offset += frame_unitcell_vectors[0]*np.round((delta[0]-offset[0])/frame_unitcell_vectors[0,0])
        for atom in anchor_molecules[next_anchor]:
            frame_positions[atom] -= offset

        # Transfer it from the available list to the used list.
        used_anchors = np.append(used_anchors, [next_anchor])
        available_anchors = np.delete(available_anchors, next_index)

    # Find the center of all anchor molecules.
    anchor_atom_indices = np.concatenate(anchor_molecules)
    center = np.mean(frame_positions[anchor_atom_indices], axis=0)
    wrap_mols(frame_positions, frame_unitcell_vectors, center, other_molecules)

def image_molecules(xyz, box, anchor_molecules, other_molecules, sorted_bonds):
    cdef vector[int[:]] omol = other_molecules
    amol = anchor_molecules
    cdef int i
    for i in range(xyz.shape[0]):
        image_frame(xyz[i], box[i], amol, omol, sorted_bonds)

def whole_molecules(xyz, box, sorted_bonds):
    cdef int i
    for i in range(xyz.shape[0]):
        make_whole(xyz[i], box[i], sorted_bonds)
