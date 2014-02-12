##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
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


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import numpy as np
from mdtraj.utils import ensure_type
import mdtraj as md
import itertools

def compute_contact_distances(traj, contacts='all', scheme='closest-heavy'):
    """
    Compute the contacts between residues in a trajectory. 

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


    Returns
    -------
    contacts : np.ndarray
        numpy array containing the contacts for each frame in the trajectory
        shape = (n_frames, n_contacts)


    See Also
    --------
    mdtraj.contact.make_square : turn the result from this function 
        into a square "contact map"
    """
    if traj.topology is None:
        raise ValueError('contact calculation requires a topology')

    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz',
                      shape=(None, None, 3), warn_on_cast=False)

    if isinstance(contacts, str):
        if contacts.lower() == 'all':
            contacts = np.array([[i, j] for i in xrange(traj.n_residues) for j in xrange(i + 3, traj.n_residues, 1)])
        else:
            raise Exception('(%s) is not a valid contacts specifier' % contacts.lower())

    else:
        if len(contacts.shape) != 2 or contacts.shape[1] != 2:
            raise Exception('contacts must be a two dimensional array with pairs of contacts in the rows')

        if contacts.max() >= traj.n_residues:
            raise Exception('contacts is requesting a residue that is not in the range of traj')

        contacts = np.array(contacts, dtype=np.int)

    
    # now the bulk of the function. This will calculate atom distances and then 
    # re-work them in the required scheme to get residue distances
    scheme = scheme.lower()
    if not scheme in ['ca', 'closest', 'closest-heavy']:
        raise Exception('scheme must be one of [ca, closest, closest-heavy]')

    if scheme == 'ca':
        res_ind_to_calpha_aind = [[atom.index for atom in residue.atoms if atom.name.lower() == 'ca'][0]
                                  for residue in traj.top.residues]
        res_ind_to_calpha_aind = np.array(res_ind_to_calpha_aind)
        # ^^^ contains the atom index for each calpha atom in each residue
        atom_pairs = res_ind_to_calpha_aind[contacts]

        ainds = np.unique(atom_pairs).flatten()
        names = [traj.top.atom(i).name for i in ainds]
        print(names)
        distances = md.compute_distances(traj, atom_pairs)
        

    elif scheme in ['closest', 'closest-heavy']:
        if scheme == 'closest':
            residue_membership = [[atom.index for atom in residue.atoms]
                                  for residue in traj.topology.residues]

        elif scheme == 'closest-heavy':
            # then remove the hydrogens from the above list
            residue_membership = [[atom.index for atom in residue.atoms if not atom.element.symbol == 'H']
                                  for residue in traj.topology.residues]

        residue_lens = [len(ainds) for ainds in residue_membership]

        atom_pairs = []
        n_atom_pairs_per_residue_pair = []
        for pair in contacts:
            atom_pairs.extend(list(itertools.product(residue_membership[pair[0]], residue_membership[pair[1]])))
            n_atom_pairs_per_residue_pair.append(residue_lens[pair[0]] * residue_lens[pair[1]])

        atom_distances = md.compute_distances(traj, atom_pairs)

        # now squash the results based on residue membership
        distances = np.zeros((len(traj), len(contacts)))
        for i in xrange(len(contacts)):
            index = np.sum(n_atom_pairs_per_residue_pair[:i])
            n = n_atom_pairs_per_residue_pair[i]
            distances[:, i] = atom_distances[:, index : index + n].min(axis=1)
        print (index + n, atom_distances.shape)

    else:
        raise ValueError('This is not supposed to happen!')

    return distances
