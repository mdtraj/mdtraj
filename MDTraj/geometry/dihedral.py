# This file is part of MDTraj.
#
# Copyright 2013 Stanford University
#
# MSMBuilder is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""Pure python code to calculate dihedral angles in a trajectory
"""
##############################################################################
# Imports
##############################################################################

import numpy as np

##############################################################################
# Functions
##############################################################################


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

    Parameters
    ----------
    top : Topology
        The topology to parse
    chain_id : int
        The index of the chain to sequence

    Notes
    -----
    By default, we assume you are interested in the first chain.
    """
    atom_dict = {}
    for chain in top.chains:
        if chain.index == chain_id:
            for residue in chain.residues:
                local_dict = {}
                for atom in residue.atoms:
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
    chain_id : int
        The index of the chain to sequence.

    Notes
    -----
    In additional finding dihedral atoms, this function could be used to
    match *general* sequences of atoms and residue_id offsets.

    Examples
    --------
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
    for chain in traj.top.chains:
        if chain.index == chain_id:
            for residue in chain.residues:
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


def strip_offsets(atom_names):
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
CHI_ATOMS = ["N", "CA", "CB", "CG"]
# CHI_ATOMS_ALT = ["N", "CA", "CB", "CG1"]  # Need to incorporate this somehow!!!

_get_indices_omega = lambda traj: atom_sequence_finder(traj, OMEGA_ATOMS)
_get_indices_phi = lambda traj: atom_sequence_finder(traj, PHI_ATOMS)
_get_indices_psi = lambda traj: atom_sequence_finder(traj, PSI_ATOMS)
_get_indices_chi = lambda traj: atom_sequence_finder(traj, CHI_ATOMS)


def compute_phi(traj):
    """Calculate the phi torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.

    Returns
    -------
    rid : np.ndarray, shape=(n_phi, 4)
        The indices of the atoms involved in each of the
        phi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_phi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_phi(traj)
    return rid, compute_dihedrals(traj, indices)


def compute_psi(traj):
    """Calculate the psi torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.

    Returns
    -------
    rid : np.ndarray, shape=(n_psi, 4)
        The indices of the atoms involved in each of the
        psi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_psi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_psi(traj)
    return rid, compute_dihedrals(traj, indices)


def compute_chi(traj):
    """Calculate the chi torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.

    Returns
    -------
    rid : np.ndarray, shape=(n_chi, 4)
        The indices of the atoms involved in each of the
        chi dihedral angles
    angles : np.ndarray, shape=(n_frames, n_chi)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_chi(traj)
    return rid, compute_dihedrals(traj, indices)


def compute_omega(traj):
    """Calculate the omega torsions of a trajectory.

    Parameters
    ----------
    traj : Trajectory
        Trajectory for which you want dihedrals.

    Returns
    -------
    rid : np.ndarray, shape=(n_omega, 4)
        The indices of the atoms involved in each of the
        omega dihedral angles
    angles : np.ndarray, shape=(n_frames, n_omega)
        The value of the dihedral angle for each of the angles in each of
        the frames.
    """
    rid, indices = _get_indices_omega(traj)
    return rid, compute_dihedrals(traj, indices)
