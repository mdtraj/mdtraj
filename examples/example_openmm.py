"""Running a simulation in OpenMM and analyzing the Ramachandran map with mdtraj.
"""

# In this example, we'e going to actually run a short simulation with OpenMM
# and analyze the results, all in a single stript.

# Obviously, running this example calculation on your machine requires
# having OpenMM installed. OpenMM can be downloaded and installed from
# https://simtk.org/home/openmm.

# Lets import some things we're going to need from mdtraj

import mdtraj
import mdtraj.reporters
import mdtraj.geometry

# And we need to get some packages from OpenMM

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

simulation.reporters.append(mdtraj.reporters.HDF5Reporter('ala2.h5', 1000))
simulation.step(1000000)
del simulation.reporters

# That was easy. Let's do a little analysis with mdtraj now.

traj = mdtraj.load('ala2.h5')
atoms, bonds = traj.topology.to_dataframe()
print atoms

# Because alanine dipeptide is a little nonstandard in the sense that
# it's basically dominated by the ACE and NME capping residues, we
# need to find the indicies of the atoms involved in the phi and psi
# angles somewhat manually. For standard cases, see geometry.compute_phi()
# and geometry.compute_psi() for easier solutions that don't require
# you to manually find the indices of each dihedral angle.

# In this case, we're just specifying the four atoms that together
# parameterize the phi and psi dihedral angles.

psi_indices, phi_indices = [6, 8, 14, 16], [4, 6, 8, 14]
angles = mdtraj.geometry.compute_dihedrals(traj, [phi_indices, psi_indices])

# Lets plot our dihedral angles in a scatter plot using matplotlib. What
# conformational states of Alanine dipeptide did we sample?

import matplotlib; matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as pp
from math import pi
pp.figure();
pp.title('Dihedral Map: Alanine dipeptide')
pp.scatter(angles[:, 0], angles[:, 1], marker='x', c=traj.time)
cbar = pp.colorbar()
cbar.set_label('Time [ps]')
pp.xlabel(r'$\Phi$ Angle [radians]')
pp.xlim(-pi, pi)
pp.ylabel(r'$\Psi$ Angle [radians]')
#@savefig ramachandran.png
pp.ylim(-pi, pi)

