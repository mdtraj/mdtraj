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


##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
# stdlib
import math
# ours
from mdtraj.core.topology import _topology_from_subset, Topology
from mdtraj.utils import unitcell
from mdtraj.utils.six import PY3
if PY3:
    basestring = str

import numpy as np

try:
    # openmm
    import openmm.unit as units
    import openmm as mm
    OPENMM_IMPORTED = True
except ImportError:
    # if someone tries to import all of mdtraj but doesn't
    # OpenMM installed, we don't want that to choke. It should
    # only choke if they actually try to USE the reporter
    OPENMM_IMPORTED = False

##############################################################################
# Classes
##############################################################################

class _BaseReporter(object):
    """
    Baseclass for reporters.
    """
    @property
    def backend(self):
        raise NotImplementedError('Must be implemented by the subclass')

    def __init__(self, file, reportInterval, coordinates=True, time=True,
                 cell=True, potentialEnergy=True, kineticEnergy=True,
                 temperature=True, velocities=False, atomSubset=None,
                 enforcePeriodicBox=None):
        """Create an OpenMM reporter

        Parameters
        ----------
        file : str, or HDF5Trajectory
            Either an open HDF5Trajecory object to write to, or a string
            specifying the filename of a new HDF5 file
        reportInterval : int
            The interval (in time steps) at which to write frames.
        coordinates : bool, default=True
            Whether to write the coordinates to the file.
        time : bool, default=True
            Whether to write the current time to the file.
        cell : bool, default=True
            Whether to write the current unitcell dimensions to the file.
        potentialEnergy : bool, default=True
            Whether to write the potential energy to the file.
        kineticEnergy : bool, default=True
            Whether to write the kinetic energy to the file.
        temperature : bool, default=True
            Whether to write the instantaneous temperature to the file.
        velocities : bool, default=False
            Whether to write the velocities of each atom to the file
        atomSubset : array_like, default=None
            Only write a subset of the atoms, with these (zero based) indices
            to the file. If None, *all* of the atoms will be written.

        Notes
        -----
        If you use the atomSubset option to write only a subset of the atoms
        to disk, the kineticEnergy, potentialEnergy, and temperature fields will
        not change. They will still refer to the energy and temperature of the *whole*
        system, and are not "subsetted" to only include the energy of your
        subsystem.
        """
        if isinstance(file, basestring):
            self._traj_file = self.backend(file, 'w')
        elif isinstance(file, self.backend):
            self._traj_file = file
            if not file.mode in ['w', 'a']:
                raise ValueError('file must be open in "w" or "a" mode')
        else:
            raise TypeError("I don't know how to handle %s" % file)

        self._reportInterval = int(reportInterval)
        self._is_intialized = False
        self._n_particles = None


        self._coordinates = bool(coordinates)
        self._time = bool(time)
        self._cell = bool(cell)
        self._potentialEnergy = bool(potentialEnergy)
        self._kineticEnergy = bool(kineticEnergy)
        self._temperature = bool(temperature)
        self._velocities = bool(velocities)
        self._needEnergy = potentialEnergy or kineticEnergy or temperature
        self._atomSubset = atomSubset
        self._atomSlice = None
        self._enforcePeriodicBox = enforcePeriodicBox

        if not OPENMM_IMPORTED:
            raise ImportError('OpenMM not found.')


    def _initialize(self, simulation):
        """Deferred initialization of the reporter, which happens before
        processing the first report.

        At the time that the first report is processed, we now have access
        to the simulation object, which we don't have at the point when the
        reporter is instantiated

        Parameters
        ----------
        simulation : openmm.app.Simulation
            The Simulation to generate a report for
        """
        if self._atomSubset is not None:
            if not min(self._atomSubset) >= 0:
                raise ValueError('atomSubset must be zero indexed. '
                                 'the smallest allowable value is zero')
            if not max(self._atomSubset) < simulation.system.getNumParticles():
                raise ValueError('atomSubset must be zero indexed. '
                                 'the largest value must be less than the number '
                                 'of particles in the system')
            if not all(a==int(a) for a in self._atomSubset):
                raise ValueError('all of the indices in atomSubset must be integers')

            self._atomSlice = self._atomSubset
            if hasattr(self._traj_file, 'topology'):
                self._traj_file.topology = _topology_from_subset(
                    Topology.from_openmm(simulation.topology), self._atomSubset)
        else:
            self._atomSlice = slice(None)
            if hasattr(self._traj_file, 'topology'):
                self._traj_file.topology = Topology.from_openmm(simulation.topology)


        system = simulation.system
        if self._temperature:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0*units.dalton:
                    dof += 3
            dof -= system.getNumConstraints()
            if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                dof -= 3
            self._dof = dof

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : openmm.app.Simulation
            The Simulation to generate a report for

        Returns
        -------
        report_description : tuple
            A five element tuple.  The first element is the number of steps
            until the next report.  The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, self._coordinates, self._velocities, False, self._needEnergy, self._enforcePeriodicBox)

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
            kwargs['time'] = state.getTime()
        if self._cell:
            vectors = state.getPeriodicBoxVectors(asNumpy=True)
            vectors = vectors.value_in_unit(getattr(units, self._traj_file.distance_unit))
            a, b, c, alpha, beta, gamma = unitcell.box_vectors_to_lengths_and_angles(*vectors)
            kwargs['cell_lengths'] = np.array([a, b, c])
            kwargs['cell_angles'] = np.array([alpha, beta, gamma])
        if self._potentialEnergy:
            kwargs['potentialEnergy'] = state.getPotentialEnergy()
        if self._kineticEnergy:
            kwargs['kineticEnergy'] = state.getKineticEnergy()
        if self._temperature:
            kwargs['temperature'] = 2*state.getKineticEnergy()/(self._dof*units.MOLAR_GAS_CONSTANT_R)
        if self._velocities:
            kwargs['velocities'] = state.getVelocities(asNumpy=True)[self._atomSlice, :]

        self._traj_file.write(*args, **kwargs)
        # flush the file to disk. it might not be necessary to do this every
        # report, but this is the most proactive solution. We don't want to
        # accumulate a lot of data in memory only to find out, at the very
        # end of the run, that there wasn't enough space on disk to hold the
        # data.
        if hasattr(self._traj_file, 'flush'):
            self._traj_file.flush()

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._needEnergy:
            energy = (state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(units.kilojoules_per_mole)
            if math.isnan(energy):
                raise ValueError('Energy is NaN')
            if math.isinf(energy):
                raise ValueError('Energy is infinite')

    def __del__(self):
        self.close()

    def close(self):
        "Close the underlying trajectory file"
        self._traj_file.close()
