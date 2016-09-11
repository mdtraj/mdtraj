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
from mdtraj.geometry import _geometry

__all__ = ['shrake_rupley']

# These van der waals radii are taken from
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# which references 
# A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001 and doi:10.1021/jp8111556. 
# M. Mantina et al. (2009). "Consistent van der Waals Radii for the Whole Main Group". J. Phys. Chem. A. 113 (19): 5806–12. doi:10.1021/jp8111556
# They also now match up with the values in current GROMACS 5.1.4
# found in share/top/vdwradii.dat.

_ATOMIC_RADII = {'H'  : 0.120, 'He' : 0.140, 'Li' : 0.182, 'Be' : 0.153,
                 'B'  : 0.192, 'C'  : 0.170, 'N'  : 0.155, 'O'  : 0.152,
                 'F'  : 0.147, 'Ne' : 0.154, 'Na' : 0.227, 'Mg' : 0.173,
                 'Al' : 0.184, 'Si' : 0.210, 'P'  : 0.180, 'S'  : 0.180,
                 'Cl' : 0.175, 'Ar' : 0.188, 'K'  : 0.275, 'Ca' : 0.231,
                 'Sc' : 0.211, 'Ni' : 0.163, 'Cu' : 0.140, 'Zn' : 0.139,
                 'Ga' : 0.187, 'Ge' : 0.211, 'As' : 0.185, 'Se' : 0.190,
                 'Br' : 0.185, 'Kr' : 0.202, 'Rb' : 0.303, 'Sr' : 0.249,
                 'Pd' : 0.163, 'Ag' : 0.172, 'Cd' : 0.158, 'In' : 0.193,
                 'Sn' : 0.217, 'Sb' : 0.206, 'Te' : 0.206, 'I'  : 0.198,
                 'Xe' : 0.216, 'Cs' : 0.343, 'Ba' : 0.268, 'Pt' : 0.175,
                 'Au' : 0.166, 'Hg' : 0.155, 'Tl' : 0.196, 'Pb' : 0.202,
                 'Bi' : 0.207, 'Po' : 0.197, 'At' : 0.202, 'Rn' : 0.220,
                 'Fr' : 0.348, 'Ra' : 0.283, 'U'  : 0.186}

##############################################################################
# Functions
##############################################################################


def shrake_rupley(traj, probe_radius=0.14, n_sphere_points=960, mode='atom'):
    """Compute the solvent accessible surface area of each atom or residue in each simulation frame.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    probe_radius : float, optional
        The radius of the probe, in nm.
    n_sphere_points : int, optional
        The number of points representing the surface of each atom, higher
        values leads to more accuracy.
    mode : {'atom', 'residue'}
        In mode == 'atom', the extracted areas are resolved per-atom
        In mode == 'residue', this is consolidated down to the
        per-residue SASA by summing over the atoms in each residue.

    Returns
    -------
    areas : np.array, shape=(n_frames, n_features)
        The accessible surface area of each atom or residue in every frame.
        If mode == 'atom', the second dimension will index the atoms in
        the trajectory, whereas if mode == 'residue', the second
        dimension will index the residues.

    Notes
    -----
    This code implements the Shrake and Rupley algorithm, with the Golden
    Section Spiral algorithm to generate the sphere points. The basic idea
    is to great a mesh of points representing the surface of each atom
    (at a distance of the van der waals radius plus the probe
    radius from the nuclei), and then count the number of such mesh points
    that are on the molecular surface -- i.e. not within the radius of another
    atom. Assuming that the points are evenly distributed, the number of points
    is directly proportional to the accessible surface area (its just 4*pi*r^2
    time the fraction of the points that are accessible).

    There are a number of different ways to generate the points on the sphere --
    possibly the best way would be to do a little "molecular dyanmics" : put the
    points on the sphere, and then run MD where all the points repel one another
    and wait for them to get to an energy minimum. But that sounds expensive.

    This code uses the golden section spiral algorithm
    (picture at http://xsisupport.com/2012/02/25/evenly-distributing-points-on-a-sphere-with-the-golden-sectionspiral/)
    where you make this spiral that traces out the unit sphere and then put points
    down equidistant along the spiral. It's cheap, but not perfect.

    The gromacs utility g_sas uses a slightly different algorithm for generating
    points on the sphere, which is based on an icosahedral tesselation.
    roughly, the icosahedral tesselation works something like this
    http://www.ziyan.info/2008/11/sphere-tessellation-using-icosahedron.html

    References
    ----------
    .. [1] Shrake, A; Rupley, JA. (1973) J Mol Biol 79 (2): 351--71.
    """

    xyz = ensure_type(traj.xyz, dtype=np.float32, ndim=3, name='traj.xyz', shape=(None, None, 3), warn_on_cast=False)
    if mode == 'atom':
        dim1 = xyz.shape[1]
        atom_mapping = np.arange(dim1, dtype=np.int32)
    elif mode == 'residue':
        dim1 = traj.n_residues
        atom_mapping = np.array(
            [a.residue.index for a in traj.top.atoms], dtype=np.int32)
        if not np.all(np.unique(atom_mapping) ==
                      np.arange(1 + np.max(atom_mapping))):
            raise ValueError('residues must have contiguous integer indices '
                             'starting from zero')
    else:
        raise ValueError('mode must be one of "residue", "atom". "%s" supplied' %
                         mode)

    out = np.zeros((xyz.shape[0], dim1), dtype=np.float32)
    atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in traj.topology.atoms]
    radii = np.array(atom_radii, np.float32) + probe_radius

    _geometry._sasa(xyz, radii, int(n_sphere_points), atom_mapping, out)

    return out
