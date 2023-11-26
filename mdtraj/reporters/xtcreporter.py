##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2018 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Joan Francesc Gilabert
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
through time in the GROMACS XTC format
"""


from __future__ import print_function, division
from ..formats.xtc import XTCTrajectoryFile
from .basereporter import _BaseReporter
from mdtraj.utils.six import PY3
if PY3:
    basestring = str
try:
    # openmm
    import openmm.unit as units
    OPENMM_IMPORTED = True
except ImportError:
    # if someone tries to import all of mdtraj but doesn't
    # OpenMM installed, we don't want that to choke. It should
    # only choke if they actually try to USE the reporter
    OPENMM_IMPORTED = False


class XTCReporter(_BaseReporter):
    """XTCReporter stores a molecular dynamics trajectory in the GROMACS
    XTC Format

    Parameters
    ----------
    file : str, or XTCTrajectoryFile
        Either an open XTCTrajectoryFile object to write to, or a string
        specifying the filename of a new XTC file to save the trajectory to.
    reportInterval : int
        The interval (in time steps) at which to write frames.
    atomSubset : array_like, default=None
        Only write a subset of the atoms, with these (zero based) indices
        to the file. If None, *all* of the atoms will be written to disk.
    append : bool, default=False
        Whether to append the trajectory to a previously existing one

    Examples
    --------
    >>> simulation = Simulation(topology, system, integrator)
    >>> xtc_reporter = XTCReporter('traj.xtc', 100)
    >>> simulation.reporters.append(xtc_reporter)
    >>> simulation.step(10000)

    >>> traj = mdtraj.trajectory.load('traj.xtc')
    """
    @property
    def backend(self):
        return XTCTrajectoryFile

    def __init__(self, file, reportInterval, atomSubset=None, append=False):
        if append:
            if isinstance(file, basestring):
                with self.backend(file, 'r') as f:
                    contents = f.read()
            elif isinstance(file, self.backend):
                raise ValueError("Currently passing an XTCTrajectoryFile in append mode is not supported, please pass a string with the filename")
            else:
                raise TypeError("I don't know how to handle %s" % file)
        super(XTCReporter, self).__init__(file, reportInterval,
            coordinates=True, time=True, cell=True, potentialEnergy=False,
            kineticEnergy=False, temperature=False, velocities=False,
            atomSubset=atomSubset)
        if append:
            self._traj_file.write(*contents)
        if not OPENMM_IMPORTED:
            raise ImportError('OpenMM not found.')

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : openmm.app.Simulation
            The Simulation to generate a report for
        state : openmm.State
            The current state of the simulation
        """
        if not self._is_intialized:
            self._initialize(simulation)
            self._is_intialized = True

        self._checkForErrors(simulation, state)
        args = ()
        kwargs = {}
        if self._coordinates:
            coordinates = state.getPositions(asNumpy=True)[self._atomSlice]
            coordinates = coordinates.value_in_unit(getattr(units, self._traj_file.distance_unit))
            args = (coordinates,)

        if self._time:
            time_step = state.getTime()
            kwargs['time'] = time_step.value_in_unit(time_step.unit)
            kwargs['step'] = simulation.currentStep
        if self._cell:
            kwargs['box'] = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(getattr(units, self._traj_file.distance_unit))
        self._traj_file.write(*args, **kwargs)
        # flush the file to disk. it might not be necessary to do this every
        # report, but this is the most proactive solution. We don't want to
        # accumulate a lot of data in memory only to find out, at the very
        # end of the run, that there wasn't enough space on disk to hold the
        # data.
        if hasattr(self._traj_file, 'flush'):
            self._traj_file.flush()
