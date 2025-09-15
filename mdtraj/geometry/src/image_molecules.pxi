
# This file is lovingly included from _geometry.pyx

import numpy as np

cimport numpy as np
from libc.stdint cimport int32_t, int64_t
from libcpp.vector cimport vector


cdef extern from "math_patch.h" nogil:
    float roundf(float x)
    float floorf(float x)

cdef int make_whole(float[:,::1] frame_positions,
                    float[:,::1] frame_unitcell_vectors,
                    int32_t[:,:] sorted_bonds) except -1 nogil:
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
        for k in range(3):
            offset[k] = frame_unitcell_vectors[2, k]*roundf(delta[2]/frame_unitcell_vectors[2,2])
        for k in range(3):
            offset[k] += frame_unitcell_vectors[1, k]*roundf((delta[1]-offset[1])/frame_unitcell_vectors[1,1])
        for k in range(3):
            offset[k] += frame_unitcell_vectors[0, k]*roundf((delta[0]-offset[0])/frame_unitcell_vectors[0,0])
            frame_positions[atom2, k] = frame_positions[atom2, k] - offset[k]

    return 0

cdef int anchor_dists(float[:,::1] frame_positions,
                      float[:,::1] frame_unitcell_vectors,
                      int[:] anchor_molecule_indices,
                      int[:] anchor_molecule_offsets,
                      float[:,:] anchor_dist,
                      int[:,:,:] anchor_nearest_atoms,
                      int num_anchors) except -1 nogil:
    cdef int ca1, ca2, mol1, mol2, i1, i2, j1, j2
    cdef float cdist
    cdef int[:] atoms1, atoms2
    for mol1 in range(num_anchors):
        j1 = anchor_molecule_offsets[mol1]
        if mol1 > 0:
            i1 = anchor_molecule_offsets[mol1-1]
        else:
            i1 = 0
        atoms1 = anchor_molecule_indices[i1:j1]
        for mol2 in range(mol1):
            j2 = anchor_molecule_offsets[mol2]
            if mol2 > 0:
                i2 = anchor_molecule_offsets[mol2-1]
            else:
                i2 = 0
            atoms2 = anchor_molecule_indices[i2:j2]
            find_closest_contact(&frame_positions[0,0], &atoms1[0], &atoms2[0], atoms1.shape[0], atoms2.shape[0], &frame_unitcell_vectors[0,0], &ca1, &ca2, &cdist)
            anchor_dist[mol1, mol2] = cdist
            anchor_dist[mol2, mol1] = cdist
            anchor_nearest_atoms[mol1, mol2, 0] = ca1
            anchor_nearest_atoms[mol1, mol2, 1] = ca2
            anchor_nearest_atoms[mol2, mol1, 0] = ca1
            anchor_nearest_atoms[mol2, mol1, 1] = ca2

    return 0

cdef int wrap_mols(float[:,::1] frame_positions,
                   float[:,::1] frame_unitcell_vectors,
                   float[:] center,
                   int[:] other_molecule_indices,
                   int[:] other_molecule_offsets) except -1 nogil:
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
    for i in range(len(other_molecule_offsets)):
        k = other_molecule_offsets[i]
        if i > 0:
            j = other_molecule_offsets[i-1]
        else:
            j = 0
        mol = other_molecule_indices[j:k]
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
        for k in range(3):
            mol_offset[k] -= frame_unitcell_vectors[1, k]*floorf(mol_offset[1]/frame_unitcell_vectors[1,1])
        for k in range(3):
            mol_offset[k] -= frame_unitcell_vectors[0, k]*floorf(mol_offset[0]/frame_unitcell_vectors[0,0])

        for j in range(n):
            for k in range(3):
                frame_positions[mol[j], k] += mol_offset[k]-mol_center[k]

    return 0

cdef void image_frame(frame_positions,
                 frame_unitcell_vectors,
                 anchor_molecule_indices,
                 anchor_molecule_offsets,
                 other_molecule_indices,
                 other_molecule_offsets,
                 sorted_bonds):

    if sorted_bonds is not None:
        make_whole(frame_positions, frame_unitcell_vectors, sorted_bonds)

    # Compute the distance between each pair of anchor molecules in this frame.
    cdef int num_anchors = len(anchor_molecule_offsets)
    anchor_dist = np.zeros((num_anchors, num_anchors), dtype=np.float32)
    cdef int[:,:,:] anchor_nearest_atoms = np.zeros((num_anchors, num_anchors, 2), dtype=np.int32)
    anchor_dists(frame_positions, frame_unitcell_vectors, anchor_molecule_indices, anchor_molecule_offsets, anchor_dist, anchor_nearest_atoms, num_anchors)


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
        j = anchor_molecule_offsets[next_anchor]
        if next_anchor > 0:
            i = anchor_molecule_offsets[next_anchor-1]
        else:
            i = 0
        if a1 in anchor_molecule_indices[i:j]:
            a2, a1 = a1, a2
        delta = frame_positions[a2]-frame_positions[a1]
        offset = frame_unitcell_vectors[2]*np.round(delta[2]/frame_unitcell_vectors[2,2])
        offset += frame_unitcell_vectors[1]*np.round((delta[1]-offset[1])/frame_unitcell_vectors[1,1])
        offset += frame_unitcell_vectors[0]*np.round((delta[0]-offset[0])/frame_unitcell_vectors[0,0])
        for atom in anchor_molecule_indices[i:j]:
            frame_positions[atom] -= offset

        # Transfer it from the available list to the used list.
        used_anchors = np.append(used_anchors, [next_anchor])
        available_anchors = np.delete(available_anchors, next_index)

    # Find the center of all anchor molecules.
    anchor_atom_indices = anchor_molecule_indices
    center = np.mean(frame_positions[anchor_atom_indices], axis=0)
    wrap_mols(frame_positions, frame_unitcell_vectors, center, other_molecule_indices, other_molecule_offsets)

def image_molecules(xyz, box, anchor_molecules, other_molecules, sorted_bonds):
    #
    # In previous versions of the code anchor_molecules and other_molecules,
    # which are lists of numpy arrays, were passed to the lower level cython
    # routines as vectors of memoryviews (vector[int[:]]), however it seems
    # some versions (?) of cython have problems with this - see e.g.
    # https://github.com/cython/cython/issues/3085
    # So here each of these lists of arrays is converted into a 1D integer
    # array of atom indices, and a second 1D integer array of pointers into
    # it, so that anchor_molecules[i] = indices[pointers[i-1]:pointers[i]]
    # and the lower level routines have been adapted to work with these, thus
    # avoiding the need for support for vectors.
    #
    n_anchors = len(anchor_molecules)
    n_others = len(other_molecules)
    anchor_molecule_offsets = np.zeros(n_anchors, dtype=np.int32)
    other_molecule_offsets = np.zeros(n_others, dtype=np.int32)
    n_anchor_atoms = 0
    for i in range(n_anchors):
        n_anchor_atoms += len(anchor_molecules[i])
        anchor_molecule_offsets[i] = n_anchor_atoms

    n_other_atoms = 0
    for i in range(n_others):
        n_other_atoms += len(other_molecules[i])
        other_molecule_offsets[i] = n_other_atoms

    anchor_molecule_indices = np.zeros(n_anchor_atoms, dtype=np.int32)
    other_molecule_indices = np.zeros(n_other_atoms, dtype=np.int32)

    offset = 0
    for i, am in enumerate(anchor_molecules):
        anchor_molecule_indices[offset:anchor_molecule_offsets[i]] = am
        offset = anchor_molecule_offsets[i]

    offset = 0
    for i, om in enumerate(other_molecules):
        other_molecule_indices[offset:other_molecule_offsets[i]] = om
        offset = other_molecule_offsets[i]

    for i in range(xyz.shape[0]):
        image_frame(xyz[i], box[i], anchor_molecule_indices,
                    anchor_molecule_offsets, other_molecule_indices,
                    other_molecule_offsets, sorted_bonds)

def whole_molecules(xyz, box, sorted_bonds):
    cdef int i
    for i in range(xyz.shape[0]):
        make_whole(xyz[i], box[i], sorted_bonds)
