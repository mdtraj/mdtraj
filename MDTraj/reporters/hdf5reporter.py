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
"""OpenMM Reporter for saving the positions of a molecular dynamics simulation
in the HDF5 format.
"""

##############################################################################
# Imports
##############################################################################
# stdlib
import math

# ours
from mdtraj.utils import unitcell
from mdtraj.hdf5 import HDF5Trajectory

import numpy as np

try:
    # openmm
    import simtk.unit as units
    import simtk.openmm as mm
    OPENMM_IMPORTED = True
except ImportError:
    # if someone tries to import all of mdtraj but doesn't
    # OpenMM installed, we don't want that to choke. It should
    # only choke if they actually try to USE the HDF5Reporter
    OPENMM_IMPORTED = False

##############################################################################
# Classes
##############################################################################


class HDF5Reporter(object):
    """HDF5Reporter stores a molecular dynamics trajectory in the HDF5 format.
    The atomic positions, periodic box vectors, and time index are saved.

    Example
    -------
    >>> simulation = Simulation(topology, system, integrator) # doctest: +SKIP
    >>> h5_reporter = HDF5Reporter('traj.h5', 100)           # doctest: +SKIP
    >>> simulation.reporters.append(h5_reporter)              # doctest: +SKIP
    >>> simulation.step(10000)                                # doctest: +SKIP

    >>> traj = mdtraj.trajectory.load('traj.lh5')             # doctest: +SKIP
    """

    def __init__(self, file, reportInterval, coordinates=True, time=True,
                 cell=True, potentialEnergy=True, kineticEnergy=True, temperature=True,
                 velocities=False):
        """Create a HDF5Reporter.

        Parameters
        ----------
        file : str, or HDF5Trajectory
            Either an open HDF5Trajecory object to write to, or a string
            specifying the filename of a new HDF5 file
        reportInterval : int
            The interval (in time steps) at which to write frames.
        coordinates : bool
            Whether to write the coordinates to the file.
        time : bool
            Whether to write the current time to the file.
        cell : bool
            Whether to write the current unitcell dimensions to the file.
        potentialEnergy : bool
            Whether to write the potential energy to the file.
        kineticEnergy : bool
            Whether to write the kinetic energy to the file.
        temperature : bool
            Whether to write the instantaneous temperature to the file.
        """
        if isinstance(file, basestring):
            self._traj_file = HDF5Trajectory(file, 'w')
        elif isinstance(file, HDF5Trajectory):
            self._traj_file = file
            if not file.mode in ['w', 'a']:
                raise ValueError('file must be open in "w" or "a" mode')
        else:
            raise TypeError("I don't know how to handle %s" % file)

        self._reportInterval = bool(reportInterval)
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
        simulation : simtk.openmm.app.Simulation
            The Simulation to generate a report for
        """
        self._traj_file.topology = simulation.topology

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

        self._traj_file.flush()

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : simtk.openmm.app.Simulation
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
        return (steps, self._coordinates, self._velocities, False, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : simtk.openmm.app.Simulation
            The Simulation to generate a report for
        state : simtk.openmm.State
            The current state of the simulation
        """
        if not self._is_intialized:
            self._initialize(simulation)
            self._is_intialized = True

        self._checkForErrors(simulation, state)

        kwargs = {}
        if self._coordinates:
            kwargs['coordinates'] = state.getPositions(asNumpy=True).astype(np.float32)
        if self._time:
            kwargs['time'] = state.getTime()
        if self._cell:
            vectors = state.getPeriodicBoxVectors(asNumpy=True)
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
            kwargs['velocities'] = state.getVelocities(asNumpy=True)

        self._traj_file.write(**kwargs)
        # flush the file to disk. it might not be necessary to do this every
        # report, but this is the most proactive solution. We don't want to
        # accumulate a lot of data in memory only to find out, at the very
        # end of the run, that there wasn't enough space on disk to hold the
        # data.
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
        self._traj_file.close()
