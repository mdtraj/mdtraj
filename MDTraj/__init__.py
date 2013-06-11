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

Examples
--------
>>> from mdtraj import trajectory                             # doctest: +SKIP
>>> t = trajectory.load('traj.dcd', topology='structure.pdb') # doctest: +SKIP
>>> t.save('traj.binpos')                                     # doctest: +SKIP

Modules
-------
- trajectory
    This is the core module, with high level interfaces to all of the loading
    and saving, as well as the central Trajectory object.
- topology
    An object oriented molecular topology.
- io
    Convenience interface to the (pytables) HDF5 library
- pdb
    PDB interface
- binpos
    Low level BINPOS interface
- dcd
    Low level DCD interface
- xtc
    Low level XTC/TRR interface
"""

__all__ = ['binpos', 'dcd', 'io', 'pdb', 'netcdf', "test", "hdf5"
            'testing', 'trajectory', 'topology', 'reporters', 'geometry']
from .hdf5 import HDF5TrajectoryFile
from .netcdf import NetCDFTrajectoryFile
from mdtraj.dcd import DCDTrajectoryFile
from .trajectory import Trajectory  # needs to be last

