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
"""The mdtraj package contains tools for loading and saving molecular dynamics
trajectories in a variety of formats, including Gromacs XTC & TRR, CHARMM/NAMD
DCD, AMBER BINPOS, PDB, and HDF5.
"""

from mdtraj.xtc import XTCTrajectoryFile
from mdtraj.trr import TRRTrajectoryFile
from .hdf5 import HDF5TrajectoryFile
from .netcdf import NetCDFTrajectoryFile
from mdtraj.dcd import DCDTrajectoryFile
from mdtraj.binpos import BINPOSTrajectoryFile
from mdtraj.pdb import PDBTrajectoryFile

from .topology import Topology
from .trajectory import *

def test():
    import nose
    nose.run()
