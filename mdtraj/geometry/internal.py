##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Kyle A. Beauchamp
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

"""Methods to calculate internal coordinates from the cartesian coordinates

This code is new and should be considered __unstable__
"""

import logging
from itertools import combinations

import numpy as np

from mdtraj.geometry.angle import compute_angles
from mdtraj.geometry.dihedral import compute_dihedrals
from mdtraj.geometry.distance import compute_distances
from mdtraj.utils import import_

# these are covalent radii taken from the crystalographic data in nm
# Dalton Trans., 2008, 2832-2838, DOI: 10.1039/B801115J
# http://pubs.rsc.org/en/Content/ArticleLanding/2008/DT/b801115j
COVALENT_RADII = {
    "C": 0.0762,
    "N": 0.0706,
    "O": 0.0661,
    "H": 0.031,
    "S": 0.105,
}
logger = logging.getLogger(__name__)

__all__ = [
    "get_redundant_internal_coordinates",
    "get_nonredundant_internal_coordinates",
    "get_connectivity",
    "get_bond_connectivity",
    "get_angle_connectivity",
    "get_dihedral_connectivity",
    "get_wilson_B",
    "get_bond_derivs",
    "get_angle_derivs",
    "get_dihedral_derivs",
]


# Get actual coordinates
def get_redundant_internal_coordinates(trajectory, **kwargs):
    """Compute internal coordinates from the cartesian coordinates

    This extracts all of the bond lengths, bond angles and dihedral angles
    from every frame in a trajectory.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        Trajectory object containing the internal coordinates

    Other Parameters
    ----------------
    ibonds : np.ndarray, optional, shape[n_bonds, 2], dtype=int
        Each row gives the indices of two atoms involved in a bond
    iangles : np.ndarray, optional shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle
    idihedrals : np.ndarray, optional, shape[n_dihedrals, 4], dtype=int
        Each row gives the indices of the four atoms which together make a
        dihedral

    Notes
    -----
    ibonds, iangles, and idihedrals will be computed usig the first
    frame in the trajectory, if not supplied

    Returns
    -------
    internal_coords : np.ndarray, shape=[n_frames, n_bonds+n_angles+n_dihedrals]
        All of the internal coordinates collected into a big array, such that
        internal_coords[i,j] gives the jth coordinate for the ith frame.
    """

    if "ibonds" in kwargs and "iangles" in kwargs and "idihedrals" in kwargs:
        ibonds = kwargs["ibonds"]
        iangles = kwargs["iangles"]
        idihedrals = kwargs["idihedrals"]
    else:
        ibonds, iangles, idihedrals = get_connectivity(trajectory)

    # convert everything to the right shape and C ordering, since
    # all of these methods are in C and are going to need things to be
    # the right type. The methods will all do a copy for things that
    # aren't the right type, but hopefully we can only do the copy once
    # instead of three times if xyzlist really does need to be reordered
    # in memory

    xyzlist = np.array(trajectory.xyz, dtype=np.float32, order="c")
    ibonds = np.array(ibonds, dtype=np.int32, order="c")
    iangles = np.array(iangles, dtype=np.int32, order="c")
    idihedrals = np.array(idihedrals, dtype=np.int32, order="c")

    b = compute_distances(xyzlist, ibonds)
    a = compute_angles(xyzlist, iangles)
    d = compute_dihedrals(xyzlist, idihedrals, degrees=False)

    return np.hstack((b, a, d))


def get_nonredundant_internal_coordinates(trajectory, conformation, get_operator=False):
    """Compute nonredudant delocalized internal coordinates from the
    cartesian coordinates

    These are basically a set of 3N-6 linear combinations of bond lengths,
    bond angles and dihedral angles that span the full space of internal
    coordinates without being redundant. The procedure to generate them
    involves collecting a bunch of "primative" internal coordinates and then
    and then taking linear combinations correspondong to eigenvectors with
    nonzero corresponding eigenvalues of G=B*B.T, where B is the so called
    "Wilson B matrix" which relates small displacements in cartesian space
    to small displacements in the internal coordinate space.

    Notes
    -----
    The projection operator from the redundant coordinate space into the
    active or nonredudant subspace is formed from the geometery in
    `conformation`, but is then applied unformly to all of the frames in
    trajectory.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        Trajectory object containing the cartesian coordinates of every
        frame in the dataset
    conformation : mdtraj.Trajectory
        Trajectort object containing a single frame (the first) to be used
        as the reference for defining the projection operator into the active
        space.
    get_operator : boolean
        Retreive the information necessary to define the cartesian ->
        nonredundant internal coordinates projection operator, including
        both the indices for generating the redudant internal coordinates
        and the linear operator that removes the redundant subspace.


    Returns
    -------
    internal_coordinates : np.ndarray, shape[n_frames, 3*N-6], dtype=float
        The position of each frame in the trajectory, represented in internal
        coordinates

    (if get_operator == True)

    activespace : np.ndarray, shape[n_redundant, n_nonredundant], dtype=float
        The linear projection operator
    ibonds : np.ndarray, shape=[n_bonds, 2], dtype=int
        n_bonds x 2 array of indices, where each row is the index of two
        atom who participate in a bond.
    iangles : np.ndarray, shape[n_angles, 3], dtype=int
        n_angles x 3 array of indices, where each row is the index of three
        atoms m,n,o such that n is bonded to both m and o.
    idihedrals : np.ndarray, shape[n_dihedrals, 4], dtype=int
        All sets of 4 atoms A,B,C,D such that A is bonded to B, B is bonded
        to C, and C is bonded to D


    References
    ----------
    Baker, Kessi, Delley J. Chem. Phys. 105, 192 (1996); doi: 10.1063/1.471864
    """
    import scipy.linalg

    ibonds, iangles, idihedrals = get_connectivity(conformation)

    B = get_wilson_B(
        conformation,
        ibonds=ibonds,
        iangles=iangles,
        idihedrals=idihedrals,
    )
    # reshape from (n_redundant, n_atoms, 3) to (n_redundant, n_atoms*3)
    B = B.reshape((B.shape[0], B.shape[1] * B.shape[2]))

    G = np.dot(B, B.T)
    eigenvalues, eigenvectors = scipy.linalg.eigh(G)

    # only the eigenvectors with nonzero eigenvalues
    # note: there should be 3N-6 of them
    activespace = eigenvectors[:, np.where(eigenvalues > 1e-10)[0]]

    if activespace.shape[1] != 3 * trajectory.xyz.shape[1] - 6:
        logger.error(
            "Active eigenspace is %dd, but 3*N - 6 = %d",
            activespace.shape[1],
            3 * trajectory.xyz.shape[1] - 6,
        )

    redundant = get_redundant_internal_coordinates(
        trajectory,
        ibonds=ibonds,
        iangles=iangles,
        idihedrals=idihedrals,
    )

    if get_operator:
        return np.dot(redundant, activespace), activespace, ibonds, iangles, idihedrals
    else:
        return np.dot(redundant, activespace)


# Comupte the connectivity, getting lists of atom indices which form bonds, bond
# angles and dihedrals
def get_connectivity(conf):
    """Get the indices of all the bonds/angles/dihedrals

    Parameters
    ----------
    conf : MDTraj.Trajectory
        An MDTraj trajectory, only the first frame will be used.

    Returns
    -------
    ibonds : np.ndarray, shape=[n_bonds, 2], dtype=int
        n_bonds x 2 array of indices, where each row is the index of two
        atom who participate in a bond.
    iangles : np.ndarray, shape[n_angles, 3], dtype=int
        n_angles x 3 array of indices, where each row is the index of three
        atoms m,n,o such that n is bonded to both m and o.
    idihedrals : np.ndarray, shape[n_dihedrals, 4], dtype=int
        All sets of 4 atoms A,B,C,D such that A is bonded to B, B is bonded
        to C, and C is bonded to D
    """

    ibonds = get_bond_connectivity(conf)
    iangles = get_angle_connectivity(ibonds)
    idihedrals = get_dihedral_connectivity(ibonds)

    return ibonds, iangles, idihedrals


def get_bond_connectivity(conf):
    """Get a list of all the bonds in a conformation

    Parameters
    ----------
    conf : MDTraj.Trajectory
        An MDTraj trajectory, only the first frame will be used.

    Returns
    -------
    ibonds : np.ndarray, shape=[n_bonds, 2], dtype=int
        n_bonds x 2 array of indices, where each row is the index of two
        atom who participate in a bond.

    Notes
    -----
    Regular bonds are assigned to all pairs of atoms where
    the interatomic distance is less than or equal to 1.3 times the
    sum of their respective covalent radii.

    References
    ----------
    Bakken and Helgaker, JCP Vol. 117, Num. 20 22 Nov. 2002
    http://folk.uio.no/helgaker/reprints/2002/JCP117b_GeoOpt.pdf
    """
    from scipy.spatial.distance import pdist, squareform

    xyz = conf.xyz[0, :, :]
    n_atoms = xyz.shape[0]

    elements = np.zeros(n_atoms, dtype="S1")
    atom_names = [a.name for a in conf.top.atoms()]
    for i in range(n_atoms):
        # name of the element that is atom[i]
        # take the first character of the AtomNames string,
        # after stripping off any digits

        elements[i] = atom_names[i].strip("123456789 ")[0]
        if elements[i] not in COVALENT_RADII.keys():
            raise ValueError(
                f"I don't know about this AtomName: {atom_names[i]}",
            )

    distance_mtx = squareform(pdist(xyz))
    connectivity = []

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Regular bonds are assigned to all pairs of atoms where
            # the interatomic distance is less than or equal to 1.3 times the
            # sum of their respective covalent radii.
            d = distance_mtx[i, j]
            if d < 1.3 * (COVALENT_RADII[elements[i]] + COVALENT_RADII[elements[j]]):
                connectivity.append((i, j))

    return np.array(connectivity)


def get_angle_connectivity(ibonds):
    """Given the bonds, get the indices of the atoms defining all the bond
    angles

    Parameters
    ----------
    ibonds : np.ndarray, shape=[n_bonds, 2], dtype=int
        n_bonds x 2 array of indices, where each row is the index of two
        atom who participate in a bond.

    Returns
    -------
    iangles : np.ndarray, shape[n_angles, 3], dtype=int
        n_angles x 3 array of indices, where each row is the index of three
        atoms m,n,o such that n is bonded to both m and o.
    """
    nx = import_("networkx")
    graph = nx.from_edgelist(ibonds)
    n_atoms = graph.number_of_nodes()
    iangles = []

    for i in range(n_atoms):
        for m, n in combinations(graph.neighbors(i), 2):
            # so now the there is a bond angle m-i-n
            iangles.append((m, i, n))

    return np.array(iangles)


def get_dihedral_connectivity(ibonds):
    """Given the bonds, get the indices of the atoms defining all the dihedral
    angles

    Parameters
    ----------
    ibonds : np.ndarray, shape=[n_bonds, 2], dtype=int
        n_bonds x 2 array of indices, where each row is the index of two
        atom who participate in a bond.

    Returns
    -------
    idihedrals : np.ndarray, shape[n_dihedrals, 4], dtype=int
        All sets of 4 atoms A,B,C,D such that A is bonded to B, B is bonded
        to C, and C is bonded to D
    """
    nx = import_("networkx")
    graph = nx.from_edgelist(ibonds)
    n_atoms = graph.number_of_nodes()
    idihedrals = []

    # TODO: CHECK FOR DIHEDRAL ANGLES THAT ARE 180 and recover
    # conf : msmbuilder.Trajectory
    #    An msmbuilder trajectory, only the first frame will be used. This
    #    is used purely to make the check for angle(ABC) != 180.

    for a in range(n_atoms):
        for b in graph.neighbors(a):
            for c in filter(lambda c: c not in [a, b], graph.neighbors(b)):
                for d in filter(lambda d: d not in [a, b, c], graph.neighbors(c)):
                    idihedrals.append((a, b, c, d))

    return np.array(idihedrals)


def get_wilson_B(conformation, **kwargs):
    """Calculate the Wilson B matrix, which collects the derivatives of the
    redundant internal coordinates w/r/t the cartesian coordinates.

    .. math::

        B_{ij} = \frac{\\partial q_i}{\\partial x_j}

    where :math:`q_i` are the internal coorindates and the :math:`x_j` are
    the Cartesian displacement coordinates of the atoms.

    BUT NOTE: THE RETURN VALUE IS ACTUALLY 3D

    Parameters
    ----------
    conformation : mdtraj.Trajectory
        Only the first frame is used

    Other Parameters
    ----------------
    ibonds : np.ndarray, optional shape[n_bonds, 2], dtype=int
        Each row gives the indices of two atoms involved in a bond
    iangles : np.ndarray, optional, shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle
    idihedrals : np.ndarray, optional, shape[n_dihedrals, 4], dtype=int
        Each row gives the indices of the four atoms which together make a
        dihedral

    Returns
    -------
    B : np.ndarray, shape=[n_internal_coordinates, n_atoms, 3]
        The layout here is 3 dimensional, where B[i,j,k] is the derivative
        of internal coordinate`q_i` with respect the cartesian coordinate which
        is the `k`-th dimension (xyz) of the `j`-th atom.
    """
    if "ibonds" in kwargs and "iangles" in kwargs and "idihedrals" in kwargs:
        ibonds = kwargs["ibonds"]
        iangles = kwargs["iangles"]
        idihedrals = kwargs["idihedrals"]
    else:
        ibonds, iangles, idihedrals = get_connectivity(conformation)

    xyz = conformation.xyz[0]

    bd = get_bond_derivs(xyz, ibonds)
    ad = get_angle_derivs(xyz, iangles)
    dd = get_dihedral_derivs(xyz, idihedrals)

    return np.vstack((bd, ad, dd))


def get_bond_derivs(xyz, ibonds):
    """Derivatives of the bond lengths with respect to cartesian coordinates

    Parameters
    ----------
    xyz : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the atomic positions for a single frame
    ibonds : np.ndarray, optional, shape[n_bonds, 2], dtype=int
        Each row gives the indices of two atoms involved in a bond

    Returns
    -------
    derivs : np.ndarray, shape=(n_bonds, n_atoms, 3)
        The gradient of the bond lengths w.r.t. each atomic position

    References
    ----------
    Bakken and Helgaker, JCP Vol. 117, Num. 20 22 Nov. 2002
    http://folk.uio.no/helgaker/reprints/2002/JCP117b_GeoOpt.pdf
    """

    n_atoms, n_bonds = xyz.shape[0], len(ibonds)

    derivatives = np.zeros((n_bonds, n_atoms, 3))
    for b, (m, n) in enumerate(ibonds):
        u = (xyz[m] - xyz[n]) / np.linalg.norm(xyz[m] - xyz[n])

        derivatives[b, m, :] = u
        derivatives[b, n, :] = -u

    return derivatives


def get_angle_derivs(xyz, iangles):
    """
    Derivatives of the bond angles with respect to cartesian coordinates

    Parameters
    ----------
    xyz : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the atomic positions for a single frame
    iangles : np.ndarray, optional shape[n_angles, 3], dtype=int
        Each row gives the indices of three atoms which together make an angle

    Returns
    -------
    derivs : np.ndarray, shape=(n_bonds, n_atoms, 3)
        The gradient of the bond angles w.r.t. each atomic position

    References
    ----------
    Bakken and Helgaker, JCP Vol. 117, Num. 20 22 Nov. 2002
    http://folk.uio.no/helgaker/reprints/2002/JCP117b_GeoOpt.pdf
    """

    n_atoms, n_angles = xyz.shape[0], len(iangles)

    derivatives = np.zeros((n_angles, n_atoms, 3))
    vector1 = np.array([1, -1, 1]) / np.sqrt(3)
    vector2 = np.array([-1, 1, 1]) / np.sqrt(3)

    for a, (m, o, n) in enumerate(iangles):
        u_prime = xyz[m] - xyz[o]
        u_norm = np.linalg.norm(u_prime)
        v_prime = xyz[n] - xyz[o]
        v_norm = np.linalg.norm(v_prime)
        u = u_prime / u_norm
        v = v_prime / v_norm

        if np.linalg.norm(u + v) < 1e-10 or np.linalg.norm(u - v) < 1e-10:
            # if they're parallel
            if np.linalg.norm(u + vector1) < 1e-10 or np.linalg.norm(u - vector1) < 1e-10:
                # and they're parallel o [1, -1, 1]
                w_prime = np.cross(u, vector2)
            else:
                w_prime = np.cross(u, vector1)
        else:
            w_prime = np.cross(u, v)

        w = w_prime / np.linalg.norm(w_prime)

        derivatives[a, m, :] = np.cross(u, w) / u_norm
        derivatives[a, n, :] = np.cross(w, v) / v_norm
        derivatives[a, o, :] = -np.cross(u, w) / u_norm - np.cross(w, v) / v_norm

    return derivatives


def get_dihedral_derivs(xyz, idihedrals):
    """
    Derivatives of the dihedral angles with respect to cartesian coordinates

    Parameters
    ----------
    xyz : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the atomic positions for a single frame
    idihedrals : np.ndarray, optional, shape[n_dihedrals, 4], dtype=int
        Each row gives the indices of the four atoms which together make a
        dihedral

    Returns
    -------
    derivs : np.ndarray, shape=(n_dihedrals, n_atoms, 3)
        The gradient of the dihedral angles w.r.t. each atomic position

    References
    ----------
    Bakken and Helgaker, JCP Vol. 117, Num. 20 22 Nov. 2002
    http://folk.uio.no/helgaker/reprints/2002/JCP117b_GeoOpt.pdf
    """

    n_atoms, n_dihedrals = xyz.shape[0], len(idihedrals)

    derivatives = np.zeros((n_dihedrals, n_atoms, 3))

    for d, (m, o, p, n) in enumerate(idihedrals):
        u_prime = xyz[m] - xyz[o]
        w_prime = xyz[p] - xyz[o]
        v_prime = xyz[n] - xyz[p]

        u_norm = np.linalg.norm(u_prime)
        w_norm = np.linalg.norm(w_prime)
        v_norm = np.linalg.norm(v_prime)

        u = u_prime / u_norm
        w = w_prime / w_norm
        v = v_prime / v_norm

        term1 = np.cross(u, w) / (u_norm * (1 - np.dot(u, w) ** 2))
        term2 = np.cross(v, w) / (v_norm * (1 - np.dot(v, w) ** 2))
        term3 = np.cross(u, w) * np.dot(u, w) / (w_norm * (1 - np.dot(u, w) ** 2))
        term4 = np.cross(v, w) * -np.dot(v, w) / (w_norm * (1 - np.dot(v, w) ** 2))

        derivatives[d, m, :] = term1
        derivatives[d, n, :] = -term2
        derivatives[d, o, :] = -term1 + term3 - term4
        derivatives[d, p, :] = term2 - term3 + term4

    return derivatives
