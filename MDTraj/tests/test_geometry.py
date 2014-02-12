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


import numpy as np

import mdtraj as md
from mdtraj.testing import get_fn, eq
from mdtraj import geometry

def test_contact():

    pdb = md.load(get_fn('bpti.pdb'))
    contacts = np.loadtxt(get_fn('contacts.dat')).astype(int)

    ca = geometry.compute_contact_distances(pdb, contacts, scheme='ca').flatten()
    closest = geometry.compute_contact_distances(pdb, contacts, scheme='closest').flatten()
    closest_heavy = geometry.compute_contact_distances(pdb, contacts, scheme='closest-heavy').flatten()

    ref_ca = np.loadtxt(get_fn('cc_ca.dat'))
    ref_closest = np.loadtxt(get_fn('cc_closest.dat'))
    ref_closest_heavy = np.loadtxt(get_fn('cc_closest-heavy.dat'))

    eq(ref_ca, ca)
    eq(ref_closest, closest)
    eq(ref_closest_heavy, closest_heavy)

    

