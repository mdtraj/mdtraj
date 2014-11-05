from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

n_steps = 10000000
output_frequency = n_steps / 400

forcefield = app.ForceField('tip3p.xml')

topology = app.Topology()
modeller = app.Modeller(topology, [])
modeller.addSolvent(forcefield, boxSize=[2.0 for i in range(3)])

topology = modeller.getTopology()
positions = modeller.getPositions()

system = forcefield.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9*u.nanometers)
integrator = mm.LangevinIntegrator(300 * u.kelvin, 1.0 / u.picoseconds, 2.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(1 * u.atmospheres, 300 * u.kelvin, 25))

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(300 * u.kelvin)
print('Equilibrating...')
simulation.step(10000)

simulation.reporters.append(app.DCDReporter('tip3p_300K_1ATM.dcd', output_frequency))
simulation.reporters.append(app.PDBReporter('tip3p_300K_1ATM.pdb', n_steps - 1))
simulation.reporters.append(app.StateDataReporter("tip3p_300K_1ATM.csv", output_frequency, step=True, temperature=True, density=True, potentialEnergy=True, separator=","))

print('Running Production...')
simulation.step(n_steps)
print('Done!')

import mdtraj as md

t = md.load("./tip3p_300K_1ATM.dcd", top="./tip3p_300K_1ATM.pdb")
t.save("./tip3p_300K_1ATM.xtc")
