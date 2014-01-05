"""Clustering with the :func:`md.rmsd` and :class:`scipy.cluster.hierarchy`.
"""

# In this example, we cluster our alanine dipeptide trajectory using
# the RMSD distance and Ward's algorithm.

import mdtraj as md
import numpy as np
import scipy.cluster.hierarchy


# Lets load up our trajectory. This is the trajectory that we generated in
# the "Running a simulation in OpenMM and analyzing the results with mdtraj"
# example. The first step is to build the rmsd cache, which precalculates
# some values for the rmsd computation.

traj = md.load('ala2.h5')

# Lets compute all pairwise rmsds between conformations.

distances = np.empty((traj.n_frames, traj.n_frames))
for i in range(traj.n_frames):
    distances[i] = md.rmsd(traj, traj, i)
print 'Max pairwise rmsd: %f nm' % np.max(distances)

# scipy.cluster implements the ward linkage
# algorithm (among others)

linkage = scipy.cluster.hierarchy.ward(distances)

# Lets plot the resulting dendrogram.

import matplotlib.pyplot as pp
pp.figure()
pp.title('RMSD Ward hierarchical clustering')
#@savefig dendrogram.png
graph = scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')
