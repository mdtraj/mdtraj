"""Principle components analysis (PCA) with :class:`scikit-learn`
"""

# scikits-learn is a premier machine learning library for python, with
# a very easy to use API and great documentation.

import mdtraj as md
from sklearn.decomposition import PCA

# Lets load up our trajectory. This is the trajectory that we generated in
# the "Running a simulation in OpenMM and analyzing the results with mdtraj"
# example.

traj = md.load('ala2.h5')
print traj

# Create a two component PCA model, and project our data down into this
# reduced dimensional space. Using just the cartesian coordinates as
# input to PCA, it's important to start with some kind of alignment.

pca1 = PCA(n_components=2)
traj.superpose(traj, 0)

reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))
print reduced_cartesian.shape

# Now we can plot the data on this projection.

import matplotlib.pyplot as pp
pp.figure()
pp.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=traj.time)
#@savefig pca_cartesian.png
pp.xlabel('PC1')
pp.ylabel('PC2')
pp.title('Cartesian coordinate PCA: alanine dipeptide')
cbar = pp.colorbar()
cbar.set_label('Time [ps]')

# Lets try cross-checking our result by using a different feature space
# that isn't sensitive to alignment, and instead to "featurize" our
# trajectory by computing the pairwise distance between every atom 
# in each frame, and using that as our high dimensional input space for PCA.

pca2 = PCA(n_components=2)
from itertools import combinations
# this python function gives you all unique pairs of elements from a list
atom_pairs = list(combinations(range(traj.n_atoms), 2))
pairwise_distances = md.geometry.compute_distances(traj, atom_pairs)
print pairwise_distances.shape
reduced_distances = pca2.fit_transform(pairwise_distances)

import matplotlib.pyplot as pp
pp.figure()
pp.scatter(reduced_distances[:, 0], reduced_distances[:,1], marker='x', c=traj.time)
#@savefig pca_distance.png
pp.xlabel('PC1')
pp.ylabel('PC2')
pp.title('Pairwise distance PCA: alanine dipeptide')
cbar = pp.colorbar()
cbar.set_label('Time [ps]')


