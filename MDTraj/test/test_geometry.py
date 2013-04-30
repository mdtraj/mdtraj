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
from mdtraj.geometry import rg
import numpy as np
from mdtraj.trajectory import load_hdf, load
import mdtraj.trajectory

def test_rg():
    t0 = load(get_fn('frame0.lh5'))
    Rg = rg.calculate_rg(t0)
    Rg0 = np.loadtxt(get_fn("Rg_frame0_ref.dat"))
    eq(Rg, Rg0)

"""
import msmbuilder.geometry.rg
from mdtraj import trajectory

r = trajectory.load("frame0.lh5")
xyz = r.xyz
Rg = msmbuilder.geometry.rg.calculate_rg(xyz)
np.savetxt("Rg_frame0_ref.dat", Rg)
"""
