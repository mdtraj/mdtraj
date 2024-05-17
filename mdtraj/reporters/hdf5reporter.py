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

"""OpenMM Reporter for saving the positions of a molecular dynamics simulation
in the HDF5 format.
"""

from mdtraj.formats.hdf5 import HDF5TrajectoryFile
from mdtraj.reporters.basereporter import _BaseReporter


class HDF5Reporter(_BaseReporter):
    """HDF5Reporter stores a molecular dynamics trajectory in the HDF5 format.

    This object supports saving all kinds of information from the simulation --
    more than any other trajectory format. In addition to all of the options,
    the topology of the system will also (of course) be stored in the file. All
    of the information is compressed, so the size of the file is not much
    different than DCD, despite the added flexibility.

    Parameters
    ----------
    file : str, or HDF5TrajectoryFile
        Either an open HDF5TrajecoryFile object to write to, or a string
        specifying the filename of a new HDF5 file to save the trajectory to.
    reportInterval : int
        The interval (in time steps) at which to write frames.
    coordinates : bool
        Whether to write the coordinates to the file.
    time : bool
        Whether to write the current time to the file.
    cell : bool
        Whether to write the current unit cell dimensions to the file.
    potentialEnergy : bool
        Whether to write the potential energy to the file.
    kineticEnergy : bool
        Whether to write the kinetic energy to the file.
    temperature : bool
        Whether to write the instantaneous temperature to the file.
    velocities : bool
        Whether to write the velocities to the file.
    atomSubset : array_like, default=None
        Only write a subset of the atoms, with these (zero based) indices
        to the file. If None, *all* of the atoms will be written to disk.
    enforcePeriodicBox: bool or None
        Specifies whether particle positions should be translated so the
        center of every molecule lies in the same periodic box. If None
        (the default), it will automatically decide whether to translate
        molecules based on whether the system being simulated uses periodic
        boundary conditions.

    Notes
    -----
    If you use the ``atomSubset`` option to write only a subset of the atoms
    to disk, the ``kineticEnergy``, ``potentialEnergy``, and ``temperature``
    fields will not change. They will still refer to the energy and temperature
    of the *whole* system, and are not "subsetted" to only include the energy
    of your subsystem.

    Examples
    --------
    >>> simulation = Simulation(topology, system, integrator)
    >>> h5_reporter = HDF5Reporter('traj.h5', 100)
    >>> simulation.reporters.append(h5_reporter)
    >>> simulation.step(10000)

    >>> traj = mdtraj.trajectory.load('traj.lh5')
    """

    @property
    def backend(self):
        return HDF5TrajectoryFile

    def __init__(
        self,
        file,
        reportInterval,
        coordinates=True,
        time=True,
        cell=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        velocities=False,
        atomSubset=None,
        enforcePeriodicBox=None,
    ):
        """Create a HDF5Reporter."""
        super().__init__(
            file,
            reportInterval,
            coordinates,
            time,
            cell,
            potentialEnergy,
            kineticEnergy,
            temperature,
            velocities,
            atomSubset,
            enforcePeriodicBox,
        )
