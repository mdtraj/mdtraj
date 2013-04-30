import numpy as np
#  http://stackoverflow.com/questions/12785836/how-to-calculate-cartesian-coordinates-from-dihedral-angle-in-python
 
def get_angle(a, b):
    """ASSUMES normalized a and b."""
    return np.arccos((a * b).sum(-1))

def normed_cross(a, b):
    v = np.cross(a, b)
    norms = (v ** 2.).sum(-1) ** 0.5
    v = (v.T / norms.T).T
    return v

def compute_dihedrals(traj, indices):
    v1 = normed_cross(x[:, indices[:,0]],x[:, indices[:,1]])
    v2 = normed_cross(x[:, indices[:,1]],x[:, indices[:,2]])
    v3 = normed_cross(x[:, indices[:,2]],x[:, indices[:,3]])

    v1v2 = normed_cross(v1, v2)
    v2v3 = normed_cross(v2, v3)
    return get_angle(v1v2, v2v3)
    
def construct_atom_dict(traj):
    """Create dictionary to lookup indices by atom name and residue_id."""
    atom_dict = {}
    for residue in traj.top.residues():
        local_dict = {}
        for atom in residue.atoms():
            local_dict[atom.name] = atom.index
        atom_dict[residue.index] = local_dict
    return atom_dict
    
def atom_sequence_finder(traj, atom_names, rid_offsets):
    atom_dict = construct_atom_dict(traj)
    num_residues = len([x for x in r.top.residues()])
    indices = []
    found_rid = []
    atoms_and_offsets = zip(atom_names, rid_offsets)
    for rid in xrange(num_residues):
        if all([atom_dict.has_key(rid + offset) for offset in rid_offsets]):  # Check that desired residue_IDs are in dict
            if all([atom_dict[rid + offset].has_key(atom) for atom, offset in atoms_and_offsets]):  # Check that we find all atom names in dict
                indices.append([atom_dict[rid + offset][atom] for atom, offset in atoms_and_offsets])  # Lookup desired atom indices and and add to list.
                found_rid.append(rid)

    indices = np.array(indices)
    found_rid = np.array(found_rid)

    return found_rid, indices

PHI_ATOMS = ["C","N","CA", "C"]
PHI_OFFSETS = [0, 1, 1, 1]
PSI_ATOMS = ["N","CA", "C", "N"]
PSI_OFFSETS = [0, 0, 0, 1]
OMEGA_ATOMS = ["CA","C","N","CA"]
OMEGA_OFFSETS = [0, 0, 1, 1]
CHI_ATOMS = ["N", "CA", "CB", "CG"]
CHI_ATOMS_ALT = ["N", "CA", "CB", "CG1"]
CHI_OFFSETS = [0, 0, 0, 0]

_get_indices_omega = lambda traj: atom_sequence_finder(traj, OMEGA_ATOMS, OMEGA_OFFSETS)
_get_indices_phi = lambda traj: atom_sequence_finder(traj, PHI_ATOMS, PHI_OFFSETS)
_get_indices_psi = lambda traj: atom_sequence_finder(traj, PSI_ATOMS, PSI_OFFSETS)
_get_indices_chi = lambda traj: atom_sequence_finder(traj, CHI_ATOMS, CHI_OFFSETS)

