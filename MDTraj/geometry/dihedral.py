# This file is part of MDTraj.
#
# Copyright 2013 Stanford University
#
# MDTraj is free software; you can redistribute it and/or modify
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

#
# Imports
#

from __future__ import print_function, division
import numpy as np
import mdtraj

from mdtraj.utils import ensure_type
try:
    import _geometry
    _HAVE_OPT = True
except ImportError:
    _HAVE_OPT = False

__all__ = ['compute_dihedrals', 'compute_phi', 'compute_psi', 'compute_omega',
           'atom_sequence_finder']

#
# Functions
#


def _dihedral(xyz, indices, out=None):
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

    return np.arctan2(p1, p2, out)


def compute_dihedrals(traj, indices, opt=True):
    """Compute the dihedral angles between the supplied quartets of atoms in each frame in a trajectory.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    indices : np.ndarray, shape=(num_dihedrals, 4), dtype=int
        Each row gives the indices of four atoms which together make a
        dihedral angle. The angle is between the planes spanned by the first
        three atoms and the last three atoms, a torsion around the bond
        between the middle two atoms.
    opt : bool, default=True
        Use an optimized native library to calculate distances. Using this
        library requires the python package "cffi" (c foreign function
        interface) which is installable via "easy_install cffi" or "pip
        install cffi". See https://pypi.python.org/pypi/cffi for more details.
        The optimized dihedral calculation is ~10-20x faster than the numpy
        implementation, so installing cffi is generally worth it.

    Returns
    -------
    dih : np.ndarray, shape=(num_dihedrals), dtype=float
        dih[i,j] gives the dihedral angle at traj[i] correponding to indices[j].

    """
    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3,
                      name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    quartets = ensure_type(np.asarray(indices), dtype=np.int32,
                           ndim=2, name='indices', shape=(None, 4), warn_on_cast=False)
    out = np.zeros((xyz.shape[0], quartets.shape[0]), dtype=np.float32)
    if _HAVE_OPT and opt:
        _geometry._dihedral(xyz, quartets, out)
    else:
        _dihedral(xyz, quartets, out)
    return out


def _construct_residue_dict(top, chain_id=0):
    res_dict = {}
    for chain in top.chains:
        if chain.index == chain_id:
            for residue in chain.residues:
                res_dict[residue.index] = residue.name
            break
    return res_dict


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
    # py3k criticial list(zip(, not just zip(, since we iterate multiple
    # times through it
    atoms_and_offsets = list(zip(atom_names, rid_offsets))
    for chain in traj.top.chains:
        if chain.index == chain_id:
            for residue in chain.residues:
                rid = residue.index
                # Check that desired residue_IDs are in dict
                if all([rid + offset in atom_dict for offset in rid_offsets]):
                    # Check that we find all atom names in dict
                    if all([atom in atom_dict[rid + offset] for atom, offset in atoms_and_offsets]):
                        # Lookup desired atom indices and and add to list.
                        atom_indices.append([atom_dict[rid + offset][atom]
                                            for atom, offset in atoms_and_offsets])
                        found_residue_ids.append(rid)

    atom_indices = np.array(atom_indices)
    found_residue_ids = np.array(found_residue_ids)
    #print(atom_indices, found_residue_ids)

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

CHI1_ATOMS = ["N", "CA", "CB", "CG"]
CHI1_ATOMS_ALT1 = ["N", "CA", "CB", "CG1"]
CHI1_ATOMS_ALT2 = ["N", "CA", "CB", "SG"]
CHI1_ATOMS_ALT3 = ["N", "CA", "CB", "OG"]
CHI1_ATOMS_ALT4 = ["N", "CA", "CB", "OG1"]


CHI2_ATOMS = ["CA", "CB", "CG", "CD"]
CHI2_ATOMS_ALT1 = ["CA", "CB", "CG", "CD1"]
CHI2_ATOMS_ALT2 = ["CA", "CB", "CG1", "CD1"]
CHI2_ATOMS_ALT3 = ["CA", "CB", "CG", "OD1"]
CHI2_ATOMS_ALT4 = ["CA", "CB", "CG", "ND1"]



CHI3_ATOMS = ["CB", "CG", "CD", "NE"]
CHI3_ATOMS_ALT1 = ["CB", "CG", "CD", "CE"]
CHI3_ATOMS_ALT2 = ["CB", "CG", "CD", "OE1"]
CHI3_ATOMS_ALT3 = ["CB", "CG", "SD", "CE"]


CHI4_ATOMS = ["CG", "CD", "NE", "CZ"]
CHI4_ATOMS_ALT = ["CG", "CD", "CE", "NZ"]


_get_indices_omega = lambda traj: atom_sequence_finder(traj, OMEGA_ATOMS)
_get_indices_phi = lambda traj: atom_sequence_finder(traj, PHI_ATOMS)
_get_indices_psi = lambda traj: atom_sequence_finder(traj, PSI_ATOMS)



_get_indices_chi1 = lambda traj: atom_sequence_finder(traj, CHI1_ATOMS)
_get_indices_chi1alt1 = lambda traj: atom_sequence_finder(traj, CHI1_ATOMS_ALT1)
_get_indices_chi1alt2 = lambda traj: atom_sequence_finder(traj, CHI1_ATOMS_ALT2)
_get_indices_chi1alt3 = lambda traj: atom_sequence_finder(traj, CHI1_ATOMS_ALT3)
_get_indices_chi1alt4 = lambda traj: atom_sequence_finder(traj, CHI1_ATOMS_ALT4)

_get_indices_chi2 = lambda traj: atom_sequence_finder(traj, CHI2_ATOMS)
_get_indices_chi2alt1 = lambda traj: atom_sequence_finder(traj, CHI2_ATOMS_ALT1)
_get_indices_chi2alt2 = lambda traj: atom_sequence_finder(traj, CHI2_ATOMS_ALT2)
_get_indices_chi2alt3 = lambda traj: atom_sequence_finder(traj, CHI2_ATOMS_ALT3)
_get_indices_chi2alt4 = lambda traj: atom_sequence_finder(traj, CHI2_ATOMS_ALT4)

_get_indices_chi3 = lambda traj: atom_sequence_finder(traj, CHI3_ATOMS)
_get_indices_chi3alt1 = lambda traj: atom_sequence_finder(traj, CHI3_ATOMS_ALT1)
_get_indices_chi3alt2 = lambda traj: atom_sequence_finder(traj, CHI3_ATOMS_ALT2)


_get_indices_chi4 = lambda traj: atom_sequence_finder(traj, CHI4_ATOMS)
_get_indices_chi4alt1 = lambda traj: atom_sequence_finder(traj, CHI4_ATOMS_ALT)


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


def compute_allchi1(traj):
    """Calculate the chi1 torsions of a trajectory.

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
    rid1, indices1 = _get_indices_chi1(traj)
    rid2, indices2 = _get_indices_chi1alt1(traj)
    rid3, indices3 = _get_indices_chi1alt2(traj)
    rid4, indices4 = _get_indices_chi1alt3(traj)
    rid5, indices5 = _get_indices_chi1alt4(traj)

    allresc1 = np.hstack((rid1, rid2,rid3,rid4,rid5))
    
    id_sort = np.argsort(allresc1)
    allresc1 = allresc1[id_sort]

    all_chi1_list = []
    [all_chi1_list.append(compute_dihedrals(traj,x)) for x in [indices1,indices2, indices3,indices4, indices5] if x.size] 
    
    all_chi1_array = np.hstack(tuple(all_chi1_list))
    all_chi1_array_sorted = all_chi1_array[:,id_sort]
    
    return allresc1,all_chi1_array_sorted   


def compute_allchi2(traj):
    """Calculate the chi2 torsions of a trajectory.

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
    rid1, indices1 = _get_indices_chi2(traj)
    rid2, indices2 = _get_indices_chi2alt1(traj)
    rid3, indices3 = _get_indices_chi2alt2(traj)
    rid4, indices4 = _get_indices_chi2alt3(traj)
    rid5, indices5 = _get_indices_chi2alt4(traj)

    allresc2 = np.hstack((rid1, rid2,rid3,rid4,rid5))
    
    id_sort = np.argsort(allresc2)
    allresc2 = allresc2[id_sort]

    all_chi2_list = []
    [all_chi2_list.append(compute_dihedrals(traj,x)) for x in [indices1, indices2, indices3,indices4,indices5] if x.size] 
    all_chi2_array = np.hstack(tuple(all_chi2_list))
    all_chi2_array_sorted = all_chi2_array[:,id_sort]

    return allresc2,all_chi2_array_sorted 


def compute_allchi3(traj):
    """Calculate the chi3 torsions of a trajectory.

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
    rid1, indices1 = _get_indices_chi3(traj)
    rid2, indices2 = _get_indices_chi3alt1(traj)
    rid3, indices3 = _get_indices_chi3alt2(traj)
    

    allresc = np.hstack((rid1, rid2,rid3))
    
    id_sort = np.argsort(allresc)
    allresc = allresc[id_sort]

    all_chi3_list = []
    [all_chi3_list.append(compute_dihedrals(traj,x)) for x in [indices1, indices2, indices3] if x.size]    
    all_chi3_array = np.hstack(tuple(all_chi3_list))
    all_chi3_array_sorted = all_chi3_array[:,id_sort]

    return allresc,all_chi3_array_sorted 

def compute_allchi4(traj):
    """Calculate the chi4 torsions of a trajectory.

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

    rid1, indices1 = _get_indices_chi4(traj)
    rid2, indices2 = _get_indices_chi4alt1(traj)
    
    allresc = np.hstack((rid1, rid2))
    
    id_sort = np.argsort(allresc)
    allresc = allresc[id_sort]

    all_chi4_list = []
    [all_chi4_list.append(compute_dihedrals(traj,x)) for x in [indices1, indices2] if x.size]                
    all_chi4_array = np.hstack(tuple(all_chi4_list))
    all_chi4_array_sorted = all_chi4_array[:,id_sort]

    return allresc,all_chi4_array_sorted 


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


if __name__ == "__main__":
    
    print("running Module")
    traj = mdtraj.load('ubq_native.pdb')
    # rc1, c1 = compute_chi1(traj)
    # rc1alt, c1alt = compute_chi1alt1(traj)

    d = _construct_residue_dict(traj.top, chain_id=0)

    resphi1,allphi = compute_phi(traj)
    respsi1,allpsi = compute_psi(traj)


    print('################')
    print(allphi.shape) 
    print(resphi1)

    print('################')
    print(allpsi.shape) 
    print(respsi1)

    reschi1,allchi1 = compute_allchi1(traj)
    print('################')
    print(allchi1.shape) 

    reschi2,allchi2 = compute_allchi2(traj)
    # compute_chi3(traj)
    print('################')
    print(allchi2.shape) 

    reschi3,allchi3 = compute_allchi3(traj)
    # compute_chi3(traj)
    print('################')
    print(allchi3.shape) 


    reschi4,allchi4 = compute_allchi4(traj)
    # compute_chi3(traj)
    print('################')
    print(allchi4.shape)



