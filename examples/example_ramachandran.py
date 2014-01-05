"""Plotting a Ramachandran map with :class:`matplotlib`
"""

import mdtraj

# Lets load up the trajectory that we simulated in a previous example

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

import matplotlib
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
