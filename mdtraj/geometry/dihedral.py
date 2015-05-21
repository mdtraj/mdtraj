##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A Beauchamp, Ravi Ramanathan, Lee-Ping Wang
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


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import _geometry, distance

__all__ = ['compute_dihedrals', 'compute_phi', 'compute_psi', 'compute_omega',
           'compute_chi1','compute_chi2','compute_chi3','compute_chi4']

##############################################################################
# Functions
##############################################################################


def _dihedral(traj, indices, periodic, out=None):
    """Compute the dihedral angles of traj for the atom indices in indices.

    Parameters
    ----------
    xyz : np.ndarray, shape=(num_frames, num_atoms, 3), dtype=float
        The XYZ coordinates of a trajectory
    indices : np.ndarray, shape=(num_dihedrals, 4), dtype=int
        Atom indices to compute dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.

    Returns
    -------
    dih : np.ndarray, shape=(num_dihedrals), dtype=float
        dih[i,j] gives the dihedral angle at traj[i] correponding to indices[j].

    """
    ix10 = indices[:, [0, 1]]
    ix21 = indices[:, [1, 2]]
    ix32 = indices[:, [2, 3]]

    b1 = distance.compute_displacements(traj, ix10, periodic=periodic, opt=False)
    b2 = distance.compute_displacements(traj, ix21, periodic=periodic, opt=False)
    b3 = distance.compute_displacements(traj, ix32, periodic=periodic, opt=False)

    c1 = np.cross(b2, b3)
    c2 = np.cross(b1, b2)

    p1 = (b1 * c1).sum(-1)
    p1 *= (b2 * b2).sum(-1) ** 0.5
    p2 = (c1 * c2).sum(-1)

    return np.arctan2(p1, p2, out)


def compute_dihedrals(traj, indices, periodic=True, opt=True):
    """Compute the dihedral angles between the supplied quartets of atoms in each frame in a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    indices : np.ndarray, shape=(n_dihedrals, 4), dtype=int
        Each row gives the indices of four atoms which together make a
        dihedral angle. The angle is between the planes spanned by the first
        three atoms and the last three atoms, a torsion around the bond
        between the middle two atoms.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    dihedrals : np.ndarray, shape=(n_frames, n_dihedrals), dtype=float
        The output array gives, in each frame from the trajectory, each of the
        `n_dihedrals` torsion angles. The angles are measured in **radians**.

    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    quartets = ensure_type(indices, dtype=np.int32, ndim=2, name='indices', shape=(None, 4), warn_on_cast=False)
    if not np.all(np.logical_and(quartets < traj.n_atoms, quartets >= 0)):
        raise ValueError('indices must be between 0 and %d' % traj.n_atoms)

    if len(quartets) == 0:
        return np.zeros((len(xyz), 0), dtype=np.float32)

    out = np.zeros((xyz.shape[0], quartets.shape[0]), dtype=np.float32)
    if periodic is True and traj._have_unitcell:
        box = ensure_type(traj.unitcell_vectors, dtype=np.float32, ndim=3, name='unitcell_vectors', shape=(len(xyz), 3, 3))
        if opt:
            _geometry._dihedral_mic(xyz, quartets, box, out)
            return out
        else:
            _dihedral(traj, quartets, periodic, out)
            return out

    if opt:
        _geometry._dihedral(xyz, quartets, out)
    else:
        _dihedral(traj, quartets, periodic, out)
    return out


def _construct_atom_dict(topology):
    """Create dictionary to lookup indices by atom name, residue_id, and chain
    index.

    Parameters
    ----------
    topology : Topology
        The topology to parse

    Returns
    -------
    atom_dict : dict
        Tree of nested dictionaries such that
        `atom_dict[chain_index][residue_index][atom_name] = atom_index`
    """
    atom_dict = {}
    for chain in topology.chains:
        residue_dict = {}
        for residue in chain.residues:
            local_dict = {}
            for atom in residue.atoms:
                local_dict[atom.name] = atom.index
            residue_dict[residue.index] = local_dict
        atom_dict[chain.index] = residue_dict

    return atom_dict


def _atom_sequence(traj, atom_names, residue_offsets=None):
    """Find sequences of atom indices corresponding to desired atoms.

    This method can be used to find sets of atoms corresponding to specific
    dihedral angles (like phi or psi). It looks for the given pattern of atoms
    in each residue of a given chain. See the example for details.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    atom_names : np.ndarray, shape=(4), dtype='str'
        Array of atoms to in each dihedral angle.
    residue_offsets : np.ndarray, optional, shape=(4), dtype='int'
        Array of integer offsets for each atom. These are used to refer
        to atoms forward or backward in the chain relative to the current
        residue

    Notes
    -----
    In additional finding dihedral atoms, this function could be used to
    match *general* sequences of atoms and residue_id offsets.

    Examples
    --------
    Here we calculate the phi torsion angles by specifying the correct
    atom names and the residue_id offsets (e.g. forward or backward in
    chain) for each atom.

    >>> traj = mdtraj.load("native.pdb")
    >>> atom_names = ["C" ,"N" , "CA", "C"]
    >>> residue_offsets = [-1, 0, 0, 0]
    >>> found_residue_ids, indices = _atom_sequence(traj, atom_names, residue_offsets)
    """
    if residue_offsets is None:
        residue_offsets = parse_offsets(atom_names)
    atom_names = _strip_offsets(atom_names)

    atom_dict = _construct_atom_dict(traj.top)

    atom_indices = []
    found_residue_ids = []
    # py3k criticial list(zip(, not just zip(, since we iterate multiple
    # times through it
    atoms_and_offsets = list(zip(atom_names, residue_offsets))
    for chain in traj.top.chains:
        cid = chain.index
        for residue in chain.residues:
            rid = residue.index
            # Check that desired residue_IDs are in dict
            if all([rid + offset in atom_dict[cid] for offset in residue_offsets]):
                # Check that we find all atom names in dict
                if all([atom in atom_dict[cid][rid + offset] for atom, offset in atoms_and_offsets]):
                    # Lookup desired atom indices and and add to list.
                    atom_indices.append([atom_dict[cid][rid + offset][atom] for atom, offset in atoms_and_offsets])
                    found_residue_ids.append(rid)

    atom_indices = np.array(atom_indices)
    found_residue_ids = np.array(found_residue_ids)

    return found_residue_ids, atom_indices


def parse_offsets(atom_names):
    """Convert a list of atom+offset strings into lists offsets.

    Parameters
    ----------
    atom_names : list
        The names of the atoms to parse for their offsets.

    Returns
    -------
    offsets : list
        The offsets of the atoms, giving whether they refer to atoms
        in the previous residue (-1), current residue (0) or next
        residue (+1)

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


def _strip_offsets(atom_names):
    """Convert a list of atom + offset strings into lists of atoms.

    Parameters
    ----------
    atom_names : list
        The names of the atoms, whose offset prexifs you want to strip

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

CHI1_ATOMS = [["N", "CA", "CB", "CG"],
              ["N", "CA", "CB", "CG1"],
              ["N", "CA", "CB", "SG"],
              ["N", "CA", "CB", "OG"],
              ["N", "CA", "CB", "OG1"]]

CHI2_ATOMS = [["CA", "CB", "CG", "CD"],
              ["CA", "CB", "CG", "CD1"],
              ["CA", "CB", "CG1", "CD1"],
              ["CA", "CB", "CG", "OD1"],
              ["CA", "CB", "CG", "ND1"],
              ["CA", "CB", "CG", "SD"]]

CHI3_ATOMS = [["CB", "CG", "CD", "NE"],
              ["CB", "CG", "CD", "CE"],
              ["CB", "CG", "CD", "OE1"],
              ["CB", "CG", "SD", "CE"]]

CHI4_ATOMS = [["CG", "CD", "NE", "CZ"],
              ["CG", "CD", "CE", "NZ"]]

_get_indices_omega = lambda traj: _atom_sequence(traj, OMEGA_ATOMS)
_get_indices_phi = lambda traj: _atom_sequence(traj, PHI_ATOMS)
_get_indices_psi = lambda traj: _atom_sequence(traj, PSI_ATOMS)


def compute_phi(traj, periodic=True, opt=True):
    """Calculate the phi torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_phi, 4)
        The indices of the atoms involved in each of the phi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_phi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_phi(traj)
    if len(indices) == 0:
        return np.empty(shape=(0, 4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)
    return indices, compute_dihedrals(traj, indices, periodic=periodic, opt=opt)


def compute_psi(traj, periodic=True, opt=True):
    """Calculate the psi torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_psi, 4)
        The indices of the atoms involved in each of the psi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_psi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_psi(traj)
    if len(indices) == 0:
        return np.empty(shape=(0, 4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)
    return indices, compute_dihedrals(traj, indices, periodic=periodic, opt=opt)


def compute_chi1(traj, periodic=True, opt=True):
    """Calculate the chi1 torsions of a trajectory. chi1 is the first side chain torsion angle
    formed between the 4 atoms over the CA-CB axis.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_chi, 4)
        The indices of the atoms involved in each of the chi1 dihedral angles
    angles : np.ndarray, shape=(n_frames, n_chi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rids, indices = zip(*(_atom_sequence(traj, atoms) for atoms in CHI1_ATOMS))
    id_sort = np.argsort(np.concatenate(rids))
    if not any(x.size for x in indices):
        return np.empty(shape=(0, 4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)

    indices = np.vstack(x for x in indices if x.size)[id_sort]
    all_chi1 = compute_dihedrals(traj, indices, periodic=periodic, opt=opt)
    return indices, all_chi1


def compute_chi2(traj, periodic=True, opt=True):
    """Calculate the chi2 torsions of a trajectory. chi2 is the second side chain torsion angle
    formed between the corresponding 4 atoms  over the CB-CG axis.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_chi, 4)
        The indices of the atoms involved in each of the chi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_chi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rids, indices = zip(*(_atom_sequence(traj, atoms) for atoms in CHI2_ATOMS))
    id_sort = np.argsort(np.concatenate(rids))
    if not any(x.size for x in indices):
        return np.empty(shape=(0, 4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)

    indices = np.vstack(x for x in indices if x.size)[id_sort]
    all_chi1 = compute_dihedrals(traj, indices, periodic=periodic, opt=opt)
    return indices, all_chi1


def compute_chi3(traj, periodic=True, opt=True):
    """Calculate the chi3 torsions of a trajectory. chi3 is the third side chain torsion angle
    formed between the corresponding 4 atoms over the CG-CD axis
    (only the residues ARG, GLN, GLU, LYS & MET have these atoms)

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_chi, 4)
        The indices of the atoms involved in each of the chi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_chi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rids, indices = zip(*(_atom_sequence(traj, atoms) for atoms in CHI3_ATOMS))
    id_sort = np.argsort(np.concatenate(rids))
    if not any(x.size for x in indices):
        return np.empty(shape=(0, 4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)

    indices = np.vstack(x for x in indices if x.size)[id_sort]
    all_chi1 = compute_dihedrals(traj, indices, periodic=periodic, opt=opt)
    return indices, all_chi1


def compute_chi4(traj, periodic=True, opt=True):
    """Calculate the chi4 torsions of a trajectory. chi4 is the fourth side chain torsion angle
    formed between the corresponding 4 atoms over the CD-CE or CD-NE axis
    (only ARG & LYS residues have these atoms)

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_chi, 4)
        The indices of the atoms involved in each of the chi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_chi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rids, indices = zip(*(_atom_sequence(traj, atoms) for atoms in CHI4_ATOMS))
    id_sort = np.argsort(np.concatenate(rids))
    if not any(x.size for x in indices):
        return np.empty(shape=(0, 4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)

    indices = np.vstack(x for x in indices if x.size)[id_sort]
    all_chi1 = compute_dihedrals(traj, indices, periodic=periodic, opt=opt)
    return indices, all_chi1


def compute_omega(traj, periodic=True, opt=True):
    """Calculate the omega torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will treat dihedrals that cross periodic images
        using the minimum image convention.
    opt : bool, default=True
        Use an optimized native library to calculate angles.

    Returns
    -------
    indices : np.ndarray, shape=(n_omega, 4)
        The indices of the atoms involved in each of the omega dihedral angles
    angles : np.ndarray, shape=(n_frames, n_omega)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_omega(traj)
    if len(indices) == 0:
        return np.empty(shape=(0,4), dtype=np.int), np.empty(shape=(len(traj), 0), dtype=np.float32)
    return indices, compute_dihedrals(traj, indices, periodic=periodic, opt=opt)
