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

from mdtraj.netcdf import NetCDFTrajectoryFile
from mdtraj.reporters.basereporter import _BaseReporter

##############################################################################
# Imports
##############################################################################

class NetCDFReporter(_BaseReporter):
    """NetCDFReporter stores a molecular dynamics trajectory in the AMBER
    NetCDF format.

    Example
    -------
    >>> simulation = Simulation(topology, system, integrator) # doctest: +SKIP
    >>> h5_reporter = NetCDFReporter('traj.h5', 100)          # doctest: +SKIP
    >>> simulation.reporters.append(h5_reporter)              # doctest: +SKIP
    >>> simulation.step(10000)                                # doctest: +SKIP

    >>> traj = mdtraj.trajectory.load('traj.nc')              # doctest: +SKIP
    """
    @property
    def backend(self):
        return NetCDFTrajectoryFile

    def __init__(self, file, reportInterval, coordinates=True, time=True,
                 cell=True, atomSubset=None):
        """Create a NetCDFReporter.

        Parameters
        ----------
        file : str, or HDF5Trajectory
            Either an open NetCDFFile object to write to, or a string
            specifying the filename of a new HDF5 file
        reportInterval : int
            The interval (in time steps) at which to write frames.
        coordinates : bool
            Whether to write the coordinates to the file.
        time : bool
            Whether to write the current time to the file.
        cell : bool
            Whether to write the current unitcell dimensions to the file.
        atomSubset : array_like, default=None
            Only write a subset of the atoms, with these (zero based) indices
            to the file. If None, *all* of the atoms will be written.
        """
        super(NetCDFReporter, self).__init__(file, reportInterval,
            coordinates, time, cell, potentialEnergy=False, kineticEnergy=False,
            temperature=False, velocities=False, atomSubset=atomSubset)
