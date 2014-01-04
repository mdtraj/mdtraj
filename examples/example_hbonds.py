"""Baker-Hubbard hydrogen bond identification
"""
import itertools
import mdtraj as md
import mdtraj.testing

# Load up some example data. This is a little 20 frame PDB, straight
# from the RCSB
path = mdtraj.testing.get_fn('2EQQ.pdb')
t = md.load(path)
print t

# :func:`md.baker_hubbard` idenfies hydrogen bonds baced on cutoffs
# for the Donor-H...Acceptor distance and angle. The criterion employed
# is :math:`\theta > 120` and :math:`r_\text{H...Acceptor} < 2.5 A` in
# at least 10% of the trajectory. The return value is a list of the 
# indices of the atoms (donor, h, acceptor) that satisfy this criteria.

hbonds = md.baker_hubbard(t)
label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
for hbond in hbonds:
    print label(hbond)

# Let's compute the actual distances between the donors and acceptors

da_distances = md.compute_distances(t, hbonds[:, [0,2]], periodic=False)

# Plot a histogram for a few of them

import matplotlib.pyplot as pp
pp.figure();
color = itertools.cycle(['r', 'b', 'gold'])
for i in [2, 3, 4]:
    pp.hist(da_distances[:, i], color=next(color), label=label(hbonds[i]), alpha=0.5)
pp.ylabel('Freq');
pp.legend();
#@savefig hbonds.png
pp.xlabel('Donor-acceptor distance [nm]')

