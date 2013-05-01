import numpy as np
# Code taken from RMG MSMBuilder2 dihedral.c

def compute_dihedrals(traj, indices):
    """Compute the dihedral angles of traj for the atom indices in indices.
    
    Parameters
    ----------
    traj : Trajectory
    indices : np.ndarray, shape=(num_dihedrals, 4), dtype=int
        
    Returns
    -------
    dih : np.ndarray, shape=(num_dihedrals), dtype=float
        dih[i,j] gives the dihedral angle at traj[i] correponding to indices[j].

    """
    x = traj.xyz
    x0 = x[:, indices[:,0]]
    x1 = x[:, indices[:,1]]
    x2 = x[:, indices[:,2]]
    x3 = x[:, indices[:,3]]    
            
    b1 = x1 - x0
    b2 = x2 - x1
    b3 = x3 - x2

    c1 = np.cross(b2, b3)
    c2 = np.cross(b1, b2)

    p1 = (b1 * c1).sum(-1)
    p1 *= (b2 * b2).sum(-1) ** 0.5
    p2 = (c1 * c2).sum(-1)
    
    return np.arctan2(p1, p2)
    
def _construct_atom_dict(traj):
    """Create dictionary to lookup indices by atom name and residue_id."""
    atom_dict = {}
    for residue in traj.top.residues():
        local_dict = {}
        for atom in residue.atoms():
            local_dict[atom.name] = atom.index
        atom_dict[residue.index] = local_dict
    return atom_dict
    
def atom_sequence_finder(traj, atom_names, rid_offsets):
    """Find sequences of atom indices correponding to desired atoms.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    atom_names : np.ndarray, shape=(4), dtype='str'
        Array of atoms to in each dihedral angle.
    rid_offsets : np.ndarray, shape=(4), dtype='int'
        Array of integer offsets for each atom.  

    Notes
    -----
    In additional finding dihedral atoms, this function could be used to
    match *general* sequences of atoms and residue_id offsets.
    """
    atom_dict = _construct_atom_dict(traj)
    indices = []
    found_rid = []
    atoms_and_offsets = zip(atom_names, rid_offsets)
    for rid in xrange(traj.n_residues):
        if all([atom_dict.has_key(rid + offset) for offset in rid_offsets]):  # Check that desired residue_IDs are in dict
            if all([atom_dict[rid + offset].has_key(atom) for atom, offset in atoms_and_offsets]):  # Check that we find all atom names in dict
                indices.append([atom_dict[rid + offset][atom] for atom, offset in atoms_and_offsets])  # Lookup desired atom indices and and add to list.
                found_rid.append(rid)

    indices = np.array(indices)
    found_rid = np.array(found_rid)

    return found_rid, indices

PHI_ATOMS = ["C","N","CA", "C"]
PHI_OFFSETS = [-1, 0, 0, 0]
PSI_ATOMS = ["N","CA", "C", "N"]
PSI_OFFSETS = [0, 0, 0, 1]
OMEGA_ATOMS = ["CA","C","N","CA"]
OMEGA_OFFSETS = [0, 0, 1, 1]
CHI_ATOMS = ["N", "CA", "CB", "CG"]
CHI_ATOMS_ALT = ["N", "CA", "CB", "CG1"]  # Need to incorporate this somehow!!!
CHI_OFFSETS = [0, 0, 0, 0]

_get_indices_omega = lambda traj: atom_sequence_finder(traj, OMEGA_ATOMS, OMEGA_OFFSETS)
_get_indices_phi = lambda traj: atom_sequence_finder(traj, PHI_ATOMS, PHI_OFFSETS)
_get_indices_psi = lambda traj: atom_sequence_finder(traj, PSI_ATOMS, PSI_OFFSETS)
_get_indices_chi = lambda traj: atom_sequence_finder(traj, CHI_ATOMS, CHI_OFFSETS)

def calculate_phi(traj):
    """Calculate the phi torsions of a trajectory."""
    rid, indices = _get_indices_phi(traj)
    return rid, compute_dihedrals(traj, indices)

def calculate_psi(traj):
    """Calculate the psi torsions of a trajectory."""
    rid, indices = _get_indices_psi(traj)
    return rid, compute_dihedrals(traj, indices)

def calculate_chi(traj):
    """Calculate the chi torsions of a trajectory."""
    rid, indices = _get_indices_chi(traj)
    return rid, compute_dihedrals(traj, indices)

def calculate_omega(traj):
    """Calculate the omega torsions of a trajectory."""
    rid, indices = _get_indices_omega(traj)
    return rid, compute_dihedrals(traj, indices)
