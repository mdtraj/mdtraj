# Copyright 2013 mdtraj developers
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
"""OpenMM Reporter for saving the state of a molecular dynamics simulation
through time in the AMBER NetCDF format
"""

##############################################################################
# Imports
##############################################################################

from mdtraj import DCDTrajectoryFile
from mdtraj.reporters.basereporter import _BaseReporter

##############################################################################
# Imports
##############################################################################

class DCDReporter(_BaseReporter):
    """DCDReporter stores a molecular dynamics trajectory in the CHARMM / NAMD
    DCD Format

    Parameters
    ----------
    file : str, or HDF5Trajectory
        Either an open HDF5Trajecory object to write to, or a string
        specifying the filename of a new HDF5 file to save the trajectory to.
    reportInterval : int
        The interval (in time steps) at which to write frames.
    atomSubset : array_like, default=None
        Only write a subset of the atoms, with these (zero based) indices
        to the file. If None, *all* of the atoms will be written to disk.

    Examples
    --------
    >>> simulation = Simulation(topology, system, integrator) # doctest: +SKIP
    >>> dcd_reporter = DCDReporter('traj.dcd', 100)           # doctest: +SKIP
    >>> simulation.reporters.append(dcd_reporter)             # doctest: +SKIP
    >>> simulation.step(10000)                                # doctest: +SKIP

    >>> traj = mdtraj.trajectory.load('traj.dcd')              # doctest: +SKIP
    """
    @property
    def backend(self):
        return DCDTrajectoryFile

    def __init__(self, file, reportInterval, atomSubset=None):
        super(DCDReporter, self).__init__(file, reportInterval,
            coordinates=True, time=False, cell=True, potentialEnergy=False,
            kineticEnergy=False, temperature=False, velocities=False,
            atomSubset=atomSubset)
