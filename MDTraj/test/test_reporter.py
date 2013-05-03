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

import os
import tempfile
import numpy as np

try:
    from simtk.unit import *
    from simtk.openmm import *
    from simtk.openmm.app import *
    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

from mdtraj.testing import get_fn, eq

from mdtraj import topology
from mdtraj import trajectory
from mdtraj.reporters import HDF5Reporter
from mdtraj.hdf5 import HDF5Trajectory

temp = tempfile.mkstemp(suffix='.h5')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file 
    this gets automatically called by nose"""
    os.unlink(temp)

@np.testing.decorators.skipif(not HAVE_OPENMM, 'No OpenMM')
def test_reporter():
    if not HAVE_OPENMM:
        return
    
    pdb = PDBFile(get_fn('native.pdb'))
    forcefield = ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic, 
        nonbondedCutoff=1.0*nanometers, constraints=HBonds, rigidWater=True)
    integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = Platform.getPlatformByName('Reference')
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    
    reporter = HDF5Reporter(temp, 2, coordinates=True, time=True,
        cell=True, potentialEnergy=True, kineticEnergy=True, temperature=True,
        velocities=True)
    simulation.reporters.append(reporter)
    simulation.step(100)
    
    reporter.close()
    
    
    with HDF5Trajectory(temp) as f:
        got = f.read()
        yield lambda: eq(got.temperature.shape, (50,))
        yield lambda: eq(got.potentialEnergy.shape, (50,))
        yield lambda: eq(got.kineticEnergy.shape, (50,))
        yield lambda: eq(got.coordinates.shape, (50, 22, 3))
        yield lambda: eq(got.velocities.shape, (50, 22, 3))
        yield lambda: eq(got.cell_lengths, 2 * np.ones((50, 3)))
        yield lambda: eq(got.cell_angles, 90*np.ones((50, 3)))
        yield lambda: eq(got.time, 0.002*2*(1+np.arange(50)))
    
        yield lambda: topology.equal(f.topology,
                                     trajectory.load(get_fn('native.pdb')).top)
        
    
