"""Rapid calculation of the RMSD drift of a simulation.
"""

import mdtraj as md

# Find two files that are distributed with MDTraj for testing purposes -- 
# we can us them to make our plot

import mdtraj.testing
crystal_fn = mdtraj.testing.get_fn('native.pdb')
trajectory_fn = mdtraj.testing.get_fn('frame0.xtc')

# Load up the trajectories from disk

crystal = md.load(crystal_fn)
trajectory = md.load(trajectory_fn, top=crystal)  # load the xtc. the crystal structure defines the topology
print trajectory


# RMSD with exchangeable hydrogen atoms is generally not a good idea
# so let's take a look at just the heavy atoms

rmsds_to_crystal = md.rmsd(trajectory, crystal, 0)
heavy_atoms = [atom.index for atom in crystal.topology.atoms if atom.element.symbol != 'H']
heavy_rmds_to_crystal = md.rmsd(trajectory, crystal, 0, atom_indices=heavy_atoms)

# Plot the results

import matplotlib.pyplot as pp
pp.figure();
pp.plot(trajectory.time, rmsds_to_crystal, 'r', label='all atom');
pp.plot(trajectory.time, heavy_rmds_to_crystal, 'b', label='heavy atom');
pp.legend();
pp.title('RMSDs to crystal');
pp.xlabel('simulation time (ps)');
#@savefig rmsds_to_crystal.png
pp.ylabel('RMSD (nm)');
