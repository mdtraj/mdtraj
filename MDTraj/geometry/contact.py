##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Christian Schwantes
# Contributors: Robert McGibbon
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
from mdtraj.utils.six.moves import xrange
from mdtraj.pdb import element
import mdtraj as md
import itertools

__all__ = ['compute_contacts', 'squareform']

##############################################################################
# Code
##############################################################################

def compute_contacts(traj, contacts='all', scheme='closest-heavy', ignore_nonprotein=True):
    """Compute the distance between pairs of residues in a trajectory.

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    contacts : np.ndarray or 'all'
        numpy array containing (0-indexed) residues to compute the
        contacts for. (e.g. np.array([[0, 10], [0, 11]]) would compute
        the contact between residue 0 and residue 10 as well as
        the contact between residue 0 and residue 11.) [NOTE: if no
        array is passed then 'all' contacts are calculated. This means
        that the result will contain all contacts between residues
        separated by at least 3 residues.]
    scheme : {'ca', 'closest', 'closest-heavy'}
        scheme to determine the distance between two residues:
            'ca' : distance between two residues is given by the distance
                between their alpha carbons
            'closest' : distance is the closest distance between any
                two atoms in the residues
            'closest-heavy' : distance is the closest distance between
                any two non-hydrogen atoms in the residues
    ignore_nonprotein : bool
        When using `contact==all`, don't compute contacts between
        "residues" which are not protein (i.e. do not contain an alpha
        carbon).

    Returns
    -------
    distances : np.ndarray, shape=(n_frames, n_pairs), dtype=np.float32
        Distances for each residue-residue contact in each frame
        of the trajectory
    residue_pairs : np.ndarray, shape=(n_pairs, 2), dtype=int
        Each row of this return value gives the indices of the residues
        involved in the contact. This argument mirrors the `contacts` input
        parameter. When `all` is specified as input, this return value
        gives the actual residue pairs resolved from `all`. Furthermore,
        when scheme=='ca', any contact pair suplied as input corresponding
        to a residue without an alpha carbon (e.g. HOH) is ignored from the
        input contacts list, meanings that the indexing of the
        output `distances` may not match up with the indexing of the input
        `contacts`. But the indexing of `distance` *will* match up with
        the indexing of `residue_pairs`

    See Also
    --------
    mdtraj.geometry.squareform : turn the result from this function
        into a square "contact map"
    Topology.residue : Get residues from the topology by index
    """
    if traj.topology is None:
        raise ValueError('contact calculation requires a topology')

    if isinstance(contacts, str):
        if contacts.lower() != 'all':
            raise ValueError('(%s) is not a valid contacts specifier' % contacts.lower())

        residue_pairs = []
        for i in xrange(traj.n_residues):
            residue_i = traj.topology.residue(i)
            if ignore_nonprotein and not any(a for a in residue_i.atoms if a.name.lower() == 'ca'):
                continue
            for j in xrange(i+3, traj.n_residues):
                residue_j = traj.topology.residue(j)
                if ignore_nonprotein and not any(a for a in residue_j.atoms if a.name.lower() == 'ca'):
                    continue
                if residue_i.chain == residue_j.chain:
                    residue_pairs.append((i, j))

        residue_pairs = np.array(residue_pairs)
        if len(residue_pairs) == 0:
            raise ValueError('No acceptable residue pairs found')

    else:
        residue_pairs = ensure_type(contacts, dtype=np.int, ndim=2, name='contacts',
                               shape=(None, 2), warn_on_cast=False)
        if not np.all((residue_pairs >= 0) * (residue_pairs < traj.n_residues)):
            raise ValueError('contacts requests a residue that is not in the permitted range')

    # now the bulk of the function. This will calculate atom distances and then
    # re-work them in the required scheme to get residue distances
    scheme = scheme.lower()
    if scheme not in ['ca', 'closest', 'closest-heavy']:
        raise ValueError('scheme must be one of [ca, closest, closest-heavy]')

    if scheme == 'ca':
        filtered_residue_pairs = []
        atom_pairs = []

        for r0, r1 in residue_pairs:
            ca_atoms_0 = [a.index for a in traj.top.residue(r0).atoms if a.name.lower() == 'ca']
            ca_atoms_1 = [a.index for a in traj.top.residue(r1).atoms if a.name.lower() == 'ca']
            if len(ca_atoms_0) == 1 and len(ca_atoms_1) == 1:
                atom_pairs.append((ca_atoms_0[0], ca_atoms_1[0]))
                filtered_residue_pairs.append((r0, r1))
            elif len(ca_atoms_0) == 0 or len(ca_atoms_1) == 0:
                # residue does not contain a CA atom, skip it
                if contacts != 'all':
                    # if the user manually asked for this residue, and didn't use "all"
                    import warnings
                    warnings.warn('Ignoring contacts pair %d-%d. No alpha carbon.' % (r0, r1))
            else:
                raise ValueError('More than 1 alpha carbon detected in residue %d or %d' % (r0, r1))

        residue_pairs = np.array(filtered_residue_pairs)
        distances = md.compute_distances(traj, atom_pairs)


    elif scheme in ['closest', 'closest-heavy']:
        if scheme == 'closest':
            residue_membership = [[atom.index for atom in residue.atoms]
                                  for residue in traj.topology.residues]
        elif scheme == 'closest-heavy':
            # then remove the hydrogens from the above list
            residue_membership = [[atom.index for atom in residue.atoms if not (atom.element == element.hydrogen)]
                                  for residue in traj.topology.residues]

        residue_lens = [len(ainds) for ainds in residue_membership]

        atom_pairs = []
        n_atom_pairs_per_residue_pair = []
        for pair in residue_pairs:
            atom_pairs.extend(list(itertools.product(residue_membership[pair[0]], residue_membership[pair[1]])))
            n_atom_pairs_per_residue_pair.append(residue_lens[pair[0]] * residue_lens[pair[1]])

        atom_distances = md.compute_distances(traj, atom_pairs)

        # now squash the results based on residue membership
        n_residue_pairs = len(residue_pairs)
        distances = np.zeros((len(traj), n_residue_pairs), dtype=np.float32)
        for i in xrange(n_residue_pairs):
            index = np.sum(n_atom_pairs_per_residue_pair[:i])
            n = n_atom_pairs_per_residue_pair[i]
            distances[:, i] = atom_distances[:, index : index + n].min(axis=1)

    else:
        raise ValueError('This is not supposed to happen!')

    return distances, residue_pairs


def squareform(distances, residue_pairs):
    """Reshape the contact distance to square contact maps

    Parameters
    ----------
    distances : np.ndarray, shape=(n_frames, n_pairs)
        Distances between pairs of residues, as computed by
        `mdtraj.geometry.compute_contacts`.
    residue_pairs : np.ndarray, shape=(n_pairs, 2)
        The indices of the residues involved in each pair, as
        returned by `mdtraj.geometry.compute_contacts`

    Returns
    -------
    contact_maps : np.ndarray, shape=(n_frames, n_residues, n_residues)
        Reshaped version of `distances`, such that the distance, in
        the `k`th frame of the trajectory from residue `i` to residue `j`
        is given by `contact_maps[k, i, j]`. All entries in `contact_maps`
        corresponding to the distance between residues that were not
        part of residue_pairs are 0.0.

    See Also
    --------
    mdtraj.compute_contacts : Compute the array of contact distances
    """
    if not isinstance(distances, np.ndarray) and distances.ndim == 2:
        raise ValueError('distances must be a 2d array')

    residue_pairs = ensure_type(
        residue_pairs, dtype=np.int, ndim=2, name='residue_pars',
        shape=(None, 2), warn_on_cast=False)

    if not np.all(residue_pairs >= 0):
        raise ValueError('residue_pairs references a residue that is not in '
                         'the permitted range')
    if distances.shape[1] != residue_pairs.shape[0]:
        raise ValueError('The number of pairs in distances, %d, does not '
                         'match the number of pairs in residue_pairs, %d.' %
                         (distances.shape[1], residue_pairs.shape[0]))

    n_residues = np.max(residue_pairs) + 1
    contact_maps = np.zeros((distances.shape[0], n_residues, n_residues), dtype=distances.dtype)
    contact_maps[:, residue_pairs[:, 0], residue_pairs[:, 1]] = distances
    contact_maps[:, residue_pairs[:, 1], residue_pairs[:, 0]] = distances

    return contact_maps
