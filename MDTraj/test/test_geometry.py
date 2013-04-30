# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

import  os
from mdtraj.testing import get_fn, eq, DocStringFormatTester
from mdtraj.geometry import rg, contact
import numpy as np
from mdtraj.trajectory import load_hdf, load
import mdtraj.trajectory

def test_rg():
    t0 = load(get_fn('frame0.lh5'))
    Rg = rg.calculate_rg(t0)
    Rg0 = np.loadtxt(get_fn("Rg_frame0_ref.dat"))
    eq(Rg, Rg0)

"""
# Compute reference data using MSMBuilder2.
import msmbuilder.geometry.rg
from mdtraj import trajectory

r = trajectory.load("frame0.lh5")
xyz = r.xyz
Rg = msmbuilder.geometry.rg.calculate_rg(xyz)
np.savetxt("Rg_frame0_ref.dat", Rg)
"""

def test_atom_distances():
    t0 = load(get_fn('frame0.lh5'))
    atom_pairs = np.loadtxt(get_fn("atom_pairs.dat"),'int')
    distances = contact.atom_distances(t0, atom_pairs)
    distances0 = np.loadtxt(get_fn("atom_distances_frame0_ref.dat"))
    eq(distances, distances0)

"""
# Compute reference data using MSMBuilder2.
import msmbuilder.geometry.contact
from mdtraj import trajectory

r = trajectory.load("frame0.lh5")
x = r.xyz
atom_pairs = np.array([[0,1],[0,2],[0,3],[1,2],[2,3],[3,4],[4,5]])
np.savetxt("./atom_pairs.dat",atom_pairs, "%d")
distances = msmbuilder.geometry.contact.atom_distances(x, atom_pairs)
np.savetxt("atom_distances_frame0_ref.dat", distances)
"""
