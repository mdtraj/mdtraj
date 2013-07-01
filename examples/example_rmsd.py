"""Plotting RMSD drift
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

# Let's also try loading only the heavy atoms from disk, because calculating
# RMSD with exchangeable hydrogen atoms is generally not a good idea

heavy_atoms = [atom.index for atom in crystal.topology.atoms if atom.element.symbol != 'H']
heavy_crystal = md.load(crystal_fn, atom_indices=heavy_atoms)
heavy_trajectory = md.load(trajectory_fn, top=crystal, atom_indices=heavy_atoms)
print heavy_trajectory

# Format the data for the RMSD calculation by building an RMSDCache.

crystal_cache = md.rmsd_cache(crystal)
trajectory_cache = md.rmsd_cache(trajectory)
crystal_heavy_cache = md.rmsd_cache(heavy_crystal)
trajectory_heavy_cache = md.rmsd_cache(heavy_trajectory)
rmsds_to_crystal = trajectory_cache.rmsds_to(crystal_cache, 0)
heavy_rmds_to_crystal = trajectory_heavy_cache.rmsds_to(crystal_heavy_cache, 0)

# Plot the results

import matplotlib.pyplot as pp
pp.figure();
pp.plot(trajectory.time, rmsds_to_crystal, 'r', label='all atom');
pp.plot(heavy_trajectory.time, heavy_rmds_to_crystal, 'b', label='heavy atom');
pp.legend();
pp.title('RMSDs to crystal');
pp.xlabel('simulation time (ps)');
#@savefig rmsds_to_crystal.png
pp.ylabel('RMSD (nm)');
