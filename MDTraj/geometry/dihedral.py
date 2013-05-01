import numpy as np


def _compute_dihedrals_xyz(xyz, indices):
    """Compute the dihedral angles of traj for the atom indices in indices.

    Parameters
    ----------
    xyz : np.ndarray, shape=(num_frames, num_atoms, 3), dtype=float
        The XYZ coordinates of a trajectory
    indices : np.ndarray, shape=(num_dihedrals, 4), dtype=int
        Atom indices to compute dihedrals.
    Returns
    -------
    dih : np.ndarray, shape=(num_dihedrals), dtype=float
        dih[i,j] gives the dihedral angle at traj[i] correponding to indices[j].

    """
    # Code adapted from RMG MSMBuilder2 dihedral.c
    x0 = xyz[:, indices[:, 0]]
    x1 = xyz[:, indices[:, 1]]
    x2 = xyz[:, indices[:, 2]]
    x3 = xyz[:, indices[:, 3]]

    b1 = x1 - x0
    b2 = x2 - x1
    b3 = x3 - x2

    c1 = np.cross(b2, b3)
    c2 = np.cross(b1, b2)

    p1 = (b1 * c1).sum(-1)
    p1 *= (b2 * b2).sum(-1) ** 0.5
    p2 = (c1 * c2).sum(-1)

    return np.arctan2(p1, p2)


def compute_dihedrals(traj, indices):
    """Compute the dihedral angles of traj for the atom indices in indices.

    Parameters
    ----------
    traj : Trajectory
        The trajectory whose dihedrals you will compute.
    indices : np.ndarray, shape=(num_dihedrals, 4), dtype=int
        Atom indices to compute dihedrals.

    Returns
    -------
    dih : np.ndarray, shape=(num_dihedrals), dtype=float
        dih[i,j] gives the dihedral angle at traj[i] correponding to indices[j].

    """
    return _compute_dihedrals_xyz(traj.xyz, indices)


def _construct_atom_dict(top, chain_id=0):
    """Create dictionary to lookup indices by atom name and residue_id.

    Notes
    -----
    By default, we assume you are interested in the first chain.
    """
    atom_dict = {}
    for chain in top.chains():
        if chain.index == chain_id:
            for residue in chain.residues():
                local_dict = {}
                for atom in residue.atoms():
                    local_dict[atom.name] = atom.index
                atom_dict[residue.index] = local_dict
            break
    return atom_dict


def atom_sequence_finder(traj, atom_names, rid_offsets=None, chain_id=0):
    """Find sequences of atom indices correponding to desired atoms.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    atom_names : np.ndarray, shape=(4), dtype='str'
        Array of atoms to in each dihedral angle.
    rid_offsets : np.ndarray, optional, shape=(4), dtype='int'
        Array of integer offsets for each atom.

    Notes
    -----
    In additional finding dihedral atoms, this function could be used to
    match *general* sequences of atoms and residue_id offsets.

    Example
    -------
    Here we calculate the phi torsion angles by specifying the correct
    atom names and the residue_id offsets (e.g. forward or backward in
    chain) for each atom.

    >>> traj = mdtraj.trajectory.load("native.pdb") # doctest: +SKIP
    >>> atom_names = ["C" ,"N" , "CA", "C"] # doctest: +SKIP
    >>> rid_offsets = [-1, 0, 0, 0] # doctest: +SKIP
    >>> found_residue_ids, indices = atom_sequence_finder(traj, atom_names, rid_offsets) # doctest: +SKIP
    """
    if rid_offsets is None:
        rid_offsets = parse_offsets(atom_names)
    atom_names = strip_offsets(atom_names)

    atom_dict = _construct_atom_dict(traj.top, chain_id=chain_id)
    atom_indices = []
    found_residue_ids = []
    atoms_and_offsets = zip(atom_names, rid_offsets)
    for chain in traj.top.chains():
        if chain.index == chain_id:
            for residue in chain.residues():
                rid = residue.index
                if all([rid + offset in atom_dict for offset in rid_offsets]):  # Check that desired residue_IDs are in dict
                    if all([atom in atom_dict[rid + offset] for atom, offset in atoms_and_offsets]):  # Check that we find all atom names in dict
                        atom_indices.append([atom_dict[rid + offset][atom] for atom, offset in atoms_and_offsets])  # Lookup desired atom indices and and add to list.
                        found_residue_ids.append(rid)

    atom_indices = np.array(atom_indices)
    found_residue_ids = np.array(found_residue_ids)

    return found_residue_ids, atom_indices

def parse_offsets(atom_names):
    """Convert a list of atom+offset strings into lists offsets.
    
    Notes
    -----
    For example, ["-C", "N", "CA", "C"] will be parsed as
    [-1, 0, 0, 0]
    """
    offsets = []
    for atom in atom_names:
        if atom[0] == "-":
            offsets.append(-1)
        elif atom[0] == "+":
            offsets.append(+1)
        else:
            offsets.append(0)
    return offsets

def strip_offsets(atom_names):
    """Convert a list of atom + offset strings into lists of atoms.
    
    Notes
    -----
    For example, ["-C", "N", "CA", "C"] will be parsed as
    ["C","N","CA","C"]
    """
    atoms = []
    for atom in atom_names:
        if atom[0] == "-":
            atoms.append(atom[1:])
        elif atom[0] == "+":
            atoms.append(atom[1:])
        else:
            atoms.append(atom)
    return atoms
    
PHI_ATOMS = ["-C", "N", "CA", "C"]
PSI_ATOMS = ["N", "CA", "C", "+N"]
OMEGA_ATOMS = ["CA", "C", "+N", "+CA"]
CHI_ATOMS = ["N", "CA", "CB", "CG"]
#CHI_ATOMS_ALT = ["N", "CA", "CB", "CG1"]  # Need to incorporate this somehow!!!

_get_indices_omega = lambda traj: atom_sequence_finder(traj, OMEGA_ATOMS)
_get_indices_phi = lambda traj: atom_sequence_finder(traj, PHI_ATOMS)
_get_indices_psi = lambda traj: atom_sequence_finder(traj, PSI_ATOMS)
_get_indices_chi = lambda traj: atom_sequence_finder(traj, CHI_ATOMS)


def compute_phi(traj):
    """Calculate the phi torsions of a trajectory."""
    rid, indices = _get_indices_phi(traj)
    return rid, compute_dihedrals(traj, indices)


def compute_psi(traj):
    """Calculate the psi torsions of a trajectory."""
    rid, indices = _get_indices_psi(traj)
    return rid, compute_dihedrals(traj, indices)


def compute_chi(traj):
    """Calculate the chi torsions of a trajectory."""
    rid, indices = _get_indices_chi(traj)
    return rid, compute_dihedrals(traj, indices)


def compute_omega(traj):
    """Calculate the omega torsions of a trajectory."""
    rid, indices = _get_indices_omega(traj)
    return rid, compute_dihedrals(traj, indices)
