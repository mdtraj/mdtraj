import numpy as np
import mdtraj.geometry.dihedral

class Featurizer(object):
    def __init__(self):
        self.keys = []
        self.n_features = 0
        self.n_feature_objects = 0
        
    def featurize(self, traj):
        return np.zeros((traj.n_frames, 0))

class VectorFeaturizer(object):
    def __init__(self, features):
        self.features = features
        self.n_features = sum([f.n_features for  f in self.features])
        self.n_feature_objects = len(f.features)

    def featurize(self, traj):
        X = []
        for f in self.features:
            X.extend(f.featurize(traj))

        return np.hstack(X)

class DihedralFeaturizer(Featurizer):
    def __init__(self, indices):
        super(Featurizer, self).__init__()
        self.indices = indices

    def featurize(self, traj):
        dih = mdtraj.geometry.dihedral.compute_dihedrals(traj, self.indices)
        return dih

class SinCosDihedralFeaturizer(Featurizer):
    def __init__(self, indices):
        super(Featurizer, self).__init__()
        self.indices = indices

    def featurize(self, traj):
        dih = mdtraj.geometry.dihedral.compute_dihedrals(traj, self.indices)
        return np.hstack((np.sin(dih), np.cos(dih)))
