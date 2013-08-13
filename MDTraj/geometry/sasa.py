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

##############################################################################
# Imports
##############################################################################

import numpy as np
from mdtraj.utils import ensure_type
from mdtraj.geometry import _HAVE_OPT
if _HAVE_OPT:
    from mdtraj.geometry import C
    from mdtraj.utils.ffi import cpointer

__all__ = ['shrake_rupley']

# these van der waals radii are taken from GROMACS 4.5.3
# and the file share/gromacs/top/vdwradii.dat
_ATOMIC_RADII = {'C': 0.15,  'F': 0.12,  'H': 0.04,
                 'N': 0.110, 'O': 0.105, 'S': 0.16}

##############################################################################
# Functions
##############################################################################


def shrake_rupley(traj, probe_radius=0.14, n_sphere_points=960):
    """Compute the solvent accessible surface area of each atom in each simulation frame.

    Parameters
    ----------
    traj : Trajectory
        An mtraj trajectory.
    probe_radius : float, optional
        The radius of the probe, in nm.
    n_sphere_pts : int, optional
        The number of points representing the surface of each atom, higher
        values leads to more accuracy.

    Returns
    -------
    areas : np.array, shape=(n_frames, n_atoms)
        The accessible surface area of each atom in every frame

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
    out = np.zeros((xyz.shape[0], xyz.shape[1]), dtype=np.float32)
    if not _HAVE_OPT:
        raise ValueError("The C library was not imported correctly")

    atom_radii = [_ATOMIC_RADII[atom.element.symbol] for atom in traj.topology.atoms]
    radii = np.array(atom_radii, np.float32) + probe_radius

    # Call kernel
    C.sasa(xyz.shape[0], xyz.shape[1], cpointer(xyz), cpointer(radii), n_sphere_points, cpointer(out))

    return out



if __name__ == '__main__':
    import mdtraj as md
    import mdtraj.testing
    t = md.load(mdtraj.testing.get_fn('2EQQ.pdb'))[0]
    out = shrake_rupley(t)
    print out
    print np.sum(out)
