"""Running a simulation in OpenMM
"""

# In this example, we'e going to actually run a short simulation with OpenMM
# saving the results to disk with MDTraj's HDF5 reporter

# Obviously, running this example calculation on your machine requires
# having OpenMM installed. OpenMM can be downloaded and installed from
# https://simtk.org/home/openmm.

# Lets import some things we're going to need from mdtraj

import os
import mdtraj
import mdtraj.reporters

# And a few things froms OpenMM

from simtk import unit
import simtk.openmm as mm
from simtk.openmm import app

# First, lets find a PDB for alanine dipeptide, the system we'll
# be simulating. We happen to have one inlcuded in the mdtraj
# package for testing named "native.pdb". Under normal circumstances
# circumstances, you shouldn't have much need for mdtraj.testing.get_fn
# (unless you're contributing tests to mdtraj!)

import mdtraj.testing
pdb = mdtraj.load(mdtraj.testing.get_fn('native.pdb'))
topology = pdb.topology.to_openmm()

# Lets use the amber99sb-ildn forcefield with implicit solvent
# and a langevin integrator. This is relatively "standard" OpenMM
# code for setting up a system.

forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic)
integrator = mm.LangevinIntegrator(330*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
simulation = app.Simulation(topology, system, integrator)

# Set the initial positions to the "first frame" of the PDB
# file (it only has one frame). Note that pdb is an mdtraj trajectory
# pass in its positions to OpenMM just fine though.

simulation.context.setPositions(pdb.xyz[0])
simulation.context.setVelocitiesToTemperature(330*unit.kelvin)

# Let's use one of the OpenMM reporters that mdtraj provides. This is
# the hdf5 reporter, which saves all kinds of information, including
# the topology, positions, energies, etc to disk. To visualize the h5
# trajectory with a non-hdf5 enabled app like PyMol or VMD, you can
# use mdconvert on the command line to easily transform it to NetCDF, DCD,
# or any other format of your preference.

if not os.path.exists('ala2.h5'):
    simulation.reporters.append(mdtraj.reporters.HDF5Reporter('ala2.h5', 1000))
    simulation.step(100000)
