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

"""OpenMM Reporter for saving the state of a molecular dynamics simulation
through time in the AMBER NetCDF format
"""

from ..formats.dcd import DCDTrajectoryFile
from .basereporter import _BaseReporter


class DCDReporter(_BaseReporter):
    """DCDReporter stores a molecular dynamics trajectory in the CHARMM / NAMD
    DCD Format

    Parameters
    ----------
    file : str, or DCDTrajectoryFile
        Either an open DCDTrajectoryFile object to write to, or a string
        specifying the filename of a new DCD file to save the trajectory to.
    reportInterval : int
        The interval (in time steps) at which to write frames.
    atomSubset : array_like, default=None
        Only write a subset of the atoms, with these (zero based) indices
        to the file. If None, *all* of the atoms will be written to disk.

    Examples
    --------
    >>> simulation = Simulation(topology, system, integrator)
    >>> dcd_reporter = DCDReporter('traj.dcd', 100)
    >>> simulation.reporters.append(dcd_reporter)
    >>> simulation.step(10000)

    >>> traj = mdtraj.trajectory.load('traj.dcd')
    """

    @property
    def backend(self):
        return DCDTrajectoryFile

    def __init__(self, file, reportInterval, atomSubset=None):
        super().__init__(
            file,
            reportInterval,
            coordinates=True,
            time=False,
            cell=True,
            potentialEnergy=False,
            kineticEnergy=False,
            temperature=False,
            velocities=False,
            atomSubset=atomSubset,
        )
