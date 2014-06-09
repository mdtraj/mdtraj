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
from itertools import product
import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry

__all__ = ['cone_wpn', 'cone_wpn_traj', 'baker_hubbard', 'kabsch_sander']

##############################################################################
# Functions
##############################################################################

def cone_wpn_(traj, exclude_water=True, periodic=False, distance_cutoff = 0.33, angle_const = 0.000044, angle_cutoff = 45):
    """Identify hydrogen bonds based on cutoffs for the Donor-H...Acceptor
    distance and angle according to the criterion outlined in [1].  
    As opposed to Baker-Hubbard, this is a "cone" criterion where the
    distance cutoff depends on the angle.

    The criterion employed is :math:`r_\\text{DA} < 3.3 A - 0.00044*\\delta_{DHA}*\\delta_{DHA}`,
    where :math:`r_\\text{DA}` is t.

    When donor the donor is 'O' and the acceptor is 'O', this corresponds to
    the definition established in [1]_. The donors considered by this method
    are NH and OH, and the acceptors considered are O. In the paper the only
    donor considered is OH.

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    periodic : bool, default=False
        Set to True to calculate displacements and angles across periodic box boundaries.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration.  Defaults to True 
        for consistency with baker_hubbard, but for all of our applications 
        it should be set to False
    distance_cutoff: float, default=0.33
        Maximum distance cutoff between heavy atoms in nanometers.
    angle_const: float, default = 0.000044
        Controls how the distance cutoff is reduced (i.e. made more strict)
        when the hydrogen bond is bent.
    angle_cutoff: float, default=45
        Any H-D-A angle greater than this number is not considered.

    Returns
    -------
    hbonds : np.array, shape=[n_hbonds, 3], dtype=int
        An array containing the indices atoms involved in each of the identified
        hydrogen bonds. Each row contains three integer indices, `(d_i, h_i,
        a_i)`, such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs (according to the definition above) in
        proportion greater than `freq` of the trajectory.

    Notes
    -----
    Each hydrogen bond is distinguished for the purpose of this function by the
    indices of the donor, hydrogen, and acceptor atoms. This means that, for
    example, when an ARG sidechain makes a hydrogen bond with its NH2 group,
    you might see what appear like double counting of the h-bonds, since the
    hydrogen bond formed via the H_1 and H_2 are counted separately, despite
    their "chemical indistinguishably"

    Examples
    --------
    >>> md.cone_wpn(t)
    array([[  0,  10,   8],
           [  0,  11,   7],
           [ 69,  73,  54],
           [ 76,  82,  65],
           [119, 131,  89],
           [140, 148, 265],
           [166, 177, 122],
           [181, 188, 231]])
    >>> label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
    >>> for hbond in hbonds:
    >>>     print label(hbond)
    GLU1-N -- GLU1-OE2 
    GLU1-N -- GLU1-OE1
    GLY6-N -- SER4-O
    CYS7-N -- GLY5-O
    TYR11-N -- VAL8-O
    MET12-N -- LYS20-O

    See Also
    --------
    baker_hubbard, kabsch_sander

    References
    ----------
    .. [1] Wernet, Ph., L.G.M. Pettersson, and A. Nilsson, et al.  
    "The Structure of the First Coordination Shell in Liquid Water." (2004) 
    Science 304, 995-999.
    """

    if traj.topology is None:
        raise ValueError('cone_wpn requires that traj contain topology '
                         'information')

    def get_donors(e0, e1):
        elems = set((e0, e1))
        bonditer = traj.topology.bonds
        atoms = [(b[0], b[1]) for b in bonditer if set((b[0].element.symbol, b[1].element.symbol)) == elems]

        indices = []
        for a0, a1 in atoms:
            if exclude_water and (a0.residue.name == 'HOH' or a1.residue.name == 'HOH'):
                continue
            pair = (a0.index, a1.index)
            # make sure to get the pair in the right order, so that the index
            # for e0 comes before e1
            if a0.element.symbol == e1:
                pair = pair[::-1]
            indices.append(pair)

        return indices

    nh_donors = get_donors('N', 'H')
    oh_donors = get_donors('O', 'H')
    xh_donors = np.array(nh_donors + oh_donors)

    if len(xh_donors) == 0:
        # if there are no hydrogens or protein in the trajectory, we get
        # no possible pairs and return nothing
        return np.zeros((0, 3), dtype=int)

    if not exclude_water:
        acceptors = [a.index for a in traj.topology.atoms if a.element.symbol == 'O']
    else:
        acceptors = [a.index for a in traj.topology.atoms if a.element.symbol == 'O' and a.residue.name != 'HOH']

    # This is used to compute the angles
    angle_triplets = np.array([(e[0][1], e[0][0], e[1]) for e in product(xh_donors, acceptors) if e[0][0] != e[1]])
    distance_pairs = angle_triplets[:, [0,2]]  # possible O..acceptor pairs

    angles = compute_angles(traj, angle_triplets, periodic=periodic) * 180.0 / np.pi # degrees
    distances = compute_distances(traj, distance_pairs, periodic=periodic, opt=True)
    cutoffs = distance_cutoff - angle_const * angles ** 2

    mask = np.logical_and(distances < cutoffs, angles < angle_cutoff)

    # This is actually returned
    angle_triplets2 = np.array([(e[0][0], e[0][1], e[1]) for e in product(xh_donors, acceptors) if e[0][0] != e[1]])

    return angle_triplets2, mask

def cone_wpn(traj, freq=0.1, exclude_water=True, periodic=False, distance_cutoff = 0.33, angle_const = 0.000044, angle_cutoff = 45):
    """ Returns a list of hydrogen bonds for the trajectory filtered by frequency. 
    See further documentation in cone_wpn_ .

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    periodic : bool, default=False
        Set to True to calculate displacements and angles across periodic box boundaries.
    freq : float, default=0.1
        Return only hydrogen bonds that occur in greater this fraction of the
        frames in the trajectory.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration.  Defaults to True 
        for consistency with baker_hubbard, but for all of our applications 
        it should be set to False

    Returns
    -------
    hbonds : np.array, shape=[n_hbonds, 3], dtype=int
        An array containing the indices atoms involved in each of the identified
        hydrogen bonds. Each row contains three integer indices, `(d_i, h_i,
        a_i)`, such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs (according to the definition above) in
        proportion greater than `freq` of the trajectory.
    """

    angle_triplets, mask = cone_wpn_(traj, exclude_water=exclude_water, periodic=periodic, distance_cutoff=distance_cutoff, angle_const=angle_const, angle_cutoff=angle_cutoff)

    # frequency of occurance of each hydrogen bond in the trajectory
    occurance = np.sum(mask, axis=0).astype(np.double) / traj.n_frames

    return angle_triplets[occurance > freq]

def cone_wpn_traj(traj, exclude_water=True, periodic=False, distance_cutoff = 0.33, angle_const = 0.000044, angle_cutoff = 45):
    """ Returns a list of hydrogen bonds for each frame in the trajectory.
    See further documentation in cone_wpn_ .

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    periodic : bool, default=False
        Set to True to calculate displacements and angles across periodic box boundaries.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration.  Defaults to True 
        for consistency with baker_hubbard, but for all of our applications 
        it should be set to False

    Returns
    -------
    hbonds : list, len=n_frames
        A list containing the atom indices involved in each of the identified
        hydrogen bonds at each frame. Each element in the list is an array
        where each row contains three integer indices, `(d_i, h_i, a_i)`, 
        such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs in that frame.
    """

    angle_triplets, mask = cone_wpn_(traj, exclude_water=exclude_water, periodic=periodic, distance_cutoff=distance_cutoff, angle_const=angle_const, angle_cutoff=angle_cutoff)
    return [angle_triplets[i] for i in mask]

def baker_hubbard(traj, freq=0.1, exclude_water=True, periodic=False):
    """Identify hydrogen bonds based on cutoffs for the Donor-H...Acceptor
    distance and angle.

    The criterion employed is :math:`\\theta > 120` and
    :math:`r_\\text{H...Acceptor} < 2.5 A`.

    When donor the donor is 'N' and the acceptor is 'O', this corresponds to
    the definition established in [1]_. The donors considered by this method
    are NH and OH, and the acceptors considered are O.

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    freq : float, default=0.1
        Return only hydrogen bonds that occur in greater this fraction of the
        frames in the trajectory.
    periodic : bool, default=False
        Set to True to calculate displacements and angles across periodic box boundaries.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration

    Returns
    -------
    hbonds : np.array, shape=[n_hbonds, 3], dtype=int
        An array containing the indices atoms involved in each of the identified
        hydrogen bonds. Each row contains three integer indices, `(d_i, h_i,
        a_i)`, such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs (according to the definition above) in
        proportion greater than `freq` of the trajectory.

    Notes
    -----
    Each hydrogen bond is distinguished for the purpose of this function by the
    indices of the donor, hydrogen, and acceptor atoms. This means that, for
    example, when an ARG sidechain makes a hydrogen bond with its NH2 group,
    you might see what appear like double counting of the h-bonds, since the
    hydrogen bond formed via the H_1 and H_2 are counted separately, despite
    their "chemical indistinguishably"

    Examples
    --------
    >>> md.baker_hubbard(t)
    array([[  0,  10,   8],
           [  0,  11,   7],
           [ 69,  73,  54],
           [ 76,  82,  65],
           [119, 131,  89],
           [140, 148, 265],
           [166, 177, 122],
           [181, 188, 231]])
    >>> label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
    >>> for hbond in hbonds:
    >>>     print label(hbond)
    GLU1-N -- GLU1-OE2 
    GLU1-N -- GLU1-OE1
    GLY6-N -- SER4-O
    CYS7-N -- GLY5-O
    TYR11-N -- VAL8-O
    MET12-N -- LYS20-O

    See Also
    --------
    cone_wpn, kabsch_sander

    References
    ----------
    .. [1] Baker, E. N., and R. E. Hubbard. "Hydrogen bonding in globular
        proteins." Progress in Biophysics and Molecular Biology
        44.2 (1984): 97-179.
    """
    # Cutoff criteria: these could be exposed as function arguments, or
    # modified if there are better definitions than the this one based only
    # on distances and angles
    distance_cutoff = 0.25            # nanometers
    angle_cutoff = 2.0 * np.pi / 3.0  # radians

    if traj.topology is None:
        raise ValueError('baker_hubbard requires that traj contain topology '
                         'information')

    def get_donors(e0, e1):
        elems = set((e0, e1))
        bonditer = traj.topology.bonds
        atoms = [(b[0], b[1]) for b in bonditer if set((b[0].element.symbol, b[1].element.symbol)) == elems]

        indices = []
        for a0, a1 in atoms:
            if exclude_water and (a0.residue.name == 'HOH' or a1.residue.name == 'HOH'):
                continue
            pair = (a0.index, a1.index)
            # make sure to get the pair in the right order, so that the index
            # for e0 comes before e1
            if a0.element.symbol == e1:
                pair = pair[::-1]
            indices.append(pair)

        return indices

    nh_donors = get_donors('N', 'H')
    oh_donors = get_donors('O', 'H')
    xh_donors = np.concatenate((nh_donors, oh_donors))

    if len(xh_donors) == 0:
        # if there are no hydrogens or protein in the trajectory, we get
        # no possible pairs and return nothing
        return np.zeros((0, 3), dtype=int)

    if not exclude_water:
        acceptors = [a.index for a in traj.topology.atoms if a.element.symbol == 'O']
    else:
        acceptors = [a.index for a in traj.topology.atoms if a.element.symbol == 'O' and a.residue.name != 'HOH']

    angle_triplets = np.array([(e[0][0], e[0][1], e[1]) for e in product(xh_donors, acceptors)])
    distance_pairs = angle_triplets[:, [1,2]]  # possible H..acceptor pairs

    angles = compute_angles(traj, angle_triplets)
    distances = compute_distances(traj, distance_pairs, periodic)

    mask = np.logical_and(distances < distance_cutoff, angles > angle_cutoff)
    # frequency of occurance of each hydrogen bond in the trajectory
    occurance = np.sum(mask, axis=0).astype(np.double) / traj.n_frames

    return angle_triplets[occurance > freq]


def kabsch_sander(traj):
    """Compute the Kabsch-Sander hydrogen bond energy between each pair
    of residues in every frame.

    Hydrogen bonds are defined using an electrostatic definition, assuming
    partial charges of -0.42 e and +0.20 e to the carbonyl oxygen and amide
    hydrogen respectively, their opposites assigned to the carbonyl carbon
    and amide nitrogen. A hydrogen bond is identified if E in the following
    equation is less than -0.5 kcal/mol:

    .. math::

        E = 0.42 \cdot 0.2 \cdot 33.2 kcal/(mol \cdot nm) * \\
            (1/r_{ON} + 1/r_{CH} - 1/r_{OH} - 1/r_{CN})

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.

    Returns
    -------
    matrices : list of scipy.sparse.csr_matrix
        The return value is a list of length equal to the number of frames
        in the trajectory. Each element is an n_residues x n_residues sparse
        matrix, where the existence of an entry at row `i`, column `j` with value
        `x` means that there exists a hydrogen bond between a backbone CO
        group at residue `i` with a backbone NH group at residue `j` whose
        Kabsch-Sander energy is less than -0.5 kcal/mol (the threshold for
        existence of the "bond"). The exact value of the energy is given by the
        value `x`.

    See Also
    --------
    cone_wpn, baker_hubbard

    References
    ----------
    .. [1] Kabsch W, Sander C (1983). "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features". Biopolymers 22 (12): 2577-637. dio:10.1002/bip.360221211
    """
    if traj.topology is None:
        raise ValueError('kabsch_sander requires topology')
    if not _geometry._processor_supports_sse41():
        raise RuntimeError('This CPU does not support the required instruction set (SSE4.1)')

    import scipy.sparse

    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz',
                      shape=(None, None, 3), warn_on_cast=False)

    ca_indices, nco_indices = [], []
    for residue in traj.topology.residues:
        if residue.name == 'PRO':
            ca_indices.append(-1)
        else:
            ca_indices.append([a.index for a in residue.atoms if a.name == 'CA'][0])

        n = [a.index for a in residue.atoms if a.name == 'N'][0]
        c = [a.index for a in residue.atoms if a.name == 'C'][0]
        o = [a.index for a in residue.atoms if a.name == 'O'][0]
        nco_indices.append([n, c, o])

    nco_indices = np.array(nco_indices, np.int32)
    ca_indices = np.array(ca_indices, np.int32)
    n_residues = len(ca_indices)
    hbonds = np.empty((xyz.shape[0], n_residues, 2), np.int32)
    henergies = np.empty((xyz.shape[0], n_residues, 2), np.float32)
    hbonds.fill(-1)
    henergies.fill(np.nan)

    _geometry._kabsch_sander(xyz, nco_indices, ca_indices, hbonds, henergies)

    # The C code returns its info in a pretty inconvenient format.
    # Let's change it to a list of scipy CSR matrices.

    matrices = []
    hbonds_mask = (hbonds != -1)
    for i in range(xyz.shape[0]):
        # appologies for this cryptic code -- we need to deal with the low
        # level aspects of the csr matrix format.
        hbonds_frame = hbonds[i]
        mask = hbonds_mask[i]
        henergies_frame = henergies[i]

        indptr = np.zeros(n_residues + 1, np.int32)
        indptr[1:] = np.cumsum(mask.sum(axis=1))
        indices = hbonds_frame[mask].flatten()
        data = henergies_frame[mask].flatten()

        matrices.append(scipy.sparse.csr_matrix((data, indices, indptr), shape=(n_residues, n_residues)))

    return matrices
