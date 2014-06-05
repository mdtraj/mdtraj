##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Lee-Ping Wang
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


from __future__ import print_function, division
import tempfile, os
import numpy as np
import mdtraj as md
from itertools import product
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.testing import get_fn, eq, DocStringFormatTester

def test_angle_nan():

    traj = md.load(get_fn("water216-eq.pdb"))

    # Cutoff criteria: these could be exposed as function arguments, or
    # modified if there are better definitions than the this one based only
    # on distances and angles
    distance_cutoff = 0.33            # nanometers
    angle_const = 0.000044            # nanometers / deg^2
    angle_cutoff = 45                 # degrees

    def get_donors(e0, e1):
        elems = set((e0, e1))
        bonditer = traj.topology.bonds
        atoms = [(b[0], b[1]) for b in bonditer if set((b[0].element.symbol, b[1].element.symbol)) == elems]

        indices = []
        for a0, a1 in atoms:
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

    acceptors = [a.index for a in traj.topology.atoms if a.element.symbol == 'O']

    angle_triplets2 = np.array([(e[0][0], e[0][1], e[1]) for e in product(xh_donors, acceptors) if e[0][0] != e[1]])

    angles = compute_angles(traj, angle_triplets2, periodic=False, opt=True) * 180.0 / np.pi # degrees

    for i in angles[0]:
        assert (not np.isnan(i))

if __name__ == "__main__":
    test_angle_nan()
