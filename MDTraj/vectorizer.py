class BaseEstimator(object):
    """Base class for all estimators"""

    def get_params(self):
        """Get parameters for this estimator.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        pass

    def set_params(self, **params):
        """Set the parameters of this estimator.

        Returns
        -------
        self
        """
        pass

    def serialize(self, ...):
        """Serialize this estimator.

        I'm not sure exactly what the API or behavior of this should be.
        It would be nice if it integrated with pytables, so that if you
        have a Pipeline with multiple estimators and/or large data
        structures, you can save it efficiently as a single file. After
        all, HDF5 is supposed to be hierachical, so you should be able to
        put multiple estimators together in a file.
        """
        return something...


class TransformerMixin(object):
    """"Mixin class for all transformers"""

    def fit_transform(self, X,  **fit_params):
        """Fit to data, then transform it.

        Parameters
        -----------
        X : any
            Generall X will be a numpy array of shape [n_samples, n_features], but
            it could be something else, like a trajectory

        Returns
        -------
        X_new : numpy array of shape [n_samples, n_features_new]
            Transformed array.
        """

        self.fit(X).transform(X)


##############################################################################
# Vectorizers: Take a trajecory as input, and return a [n_frames, n_features]
# array
##############################################################################


class PositionVectorizer(BaseEstimator, TransformerMixin):
    """Transforms a trajectory into a set of pairwise distances

    This transformer will extract the positions of the atoms in a trajectory
    after aliging (RMSD) the structure to a reference frame.

    Parameters
    ----------
    reference : Trajectry
        The structure with which to align the fit trajectories to. If
        `reference` is a trajectory containing more than one frame, only
        the first frame will be used to align structures to.
    alignment_indices : numpy array of shape [n_alignment_indices]
        The indices of the atoms with which to perform RMSD alignment. If
        both `reference` and `alignment_indices` are none, no alignment
        will be done, and the cartesian coordinates of the trajectories
        will be used *as is*.
    """

    def __init__(self, reference=None, alignment_indices=None):
        self.reference = reference[0]
        self.alignment_indices = alignment_indices

    def fit(self, X):
        """This is a no-op, because PositionVectorizer does not need to
        be fit to training data.

        Returns
        -------
        self
        """
        return self

    def transform(self, X):
        """Extract the positions of all of the atoms in a trajectory
        after aliging them to a reference

        Parameters
        ----------
        X : md.Trajectory
            A molecular dynamics trajectory

        Returns
        -------
        X_new : numpy array of shape [n_frames, n_atoms*3]
            `X_new[i, 3*j+d]` will contain the cartesian coordinate of the
            `i`-th frame in the `d`th dimension (x, y or z) after alignment.
        """
        # do the alignment
        return self


    def inverse_transform(self, X):


class DistanceVectorizer(BaseEstimator, TransformerMixin):
    """Transforms a trajectory into a set of pairwise distances

    This transformer turns trajectories into vectors of pairwise distances
    between specified atoms.

    Parameters
    ----------
    pair_indices : numpy_array of shape [n_distances, 2]
        Pairs of indices of the atoms between which you wish to calculate
        distances
    use_periodic_boundries : bool
        Compute distances accross periodic boundary conditions. This is used
        only when when the trajectories contain PBC information.

    Examples
    --------
    >>> X = md.load('trajectory.h5')
    >>> distances = DistanceVectorizer([[1,2]]).fit_transform(X)
    """

    def __init__(self, pair_indices, use_periodic_boundries=False):
        self.pair_indices = pair_indices
        self.use_periodic_boundries = use_periodic_boundries

    def fit(self, X):
        """This is a no-op, because DistanceVectorizer does not need to
        be fit to training data.

        Returns
        -------
        self
        """
        return self

    def transform(self, X):
        """Extract the distance between pairs of atoms from a trajectory

        Parameters
        ----------
        X : md.Trajectory
            A molecular dynamics trajectory

        Returns
        -------
        d : numpy array of shape [n_frames, n_distances]
        """
        # ... calculate the distances
        return distances


class AngleVectorizer(BaseEstimator, TransformerMixin):
    """Transforms a trajectory into a set of angles

    This transformer turns trajectories into vectors of angles
    between specified triplets atoms.

    Parameters
    ----------
    triplet_indices : numpy_array of shape [n_angles, 3]
        Each row of specified three indices, p0, p1, p2 of atoms. The
        calculated angle will be around the central atom, p1.

    Notes
    -----
    This vectorizer will not inspect periodic boundary conditions, and
    will give incorrect results if an angle bridges across periodic images.

    Examples
    --------
    >>> X = md.load('trajectory.h5')
    >>> distances = AngleVectorizer([[0, 1, 2]]).fit_transform(X)
    """

    def __init__(self, triplet_indices):
        self.triplet_indices = triplet_indices

    def fit(self, X):
        """This is a no-op, because AngleVectorizer does not need to
        be fit to training data.

        Returns
        -------
        self
        """
        return self

    def transform(self, X):
        """Extract the angles between atoms from a trajectory

        Parameters
        ----------
        X : md.Trajectory
            A molecular dynamics trajectory

        Returns
        -------
        d : numpy array of shape [n_frames, n_angles]
            The angles, in radians
        """
        # ... calculate the angles
        return angles


class DihedralVectorizer(BaseEstimator, TransformerMixin):
    """Transforms a trajectory into a set of torsions

    This transformer turns trajectories into vectors of torsion angles
    between specified quartets of atoms.

    Notes
    -----
    This vectorizer will not inspect periodic boundary conditions, and
    will give incorrect results if an angle bridges across periodic images.

    Parameters
    ----------
    quartet_indices_indices : numpy_array of shape [n_dihedrals, 4]
        Each row of specifies four indices, p0, p1, p2, p3 of atoms. The
        calculated angle will be between the plane formed by atoms p0, p1,
        and p2 and the plane formed by atoms p1, p2, and p3.

    Examples
    --------
    >>> X = md.load('trajectory.h5')
    >>> distances = DihedralVectorizer([[0, 1, 2, 3]]).fit_transform(X)
    """
    def __init__(self, quartet_indices):
        self.quartet_indices = quartet_indices

    def fit(self, X):
        """This is a no-op, because DihedralVectorizer does not need to
        be fit to training data.

        Returns
        -------
        self
        """
        return self

    def transform(X):
        """Extract torsion angles from a trajectory

        Parameters
        ----------
        X : md.Trajectory
            A molecular dynamics trajectory

        Returns
        -------
        d : numpy array of shape [n_frames, n_dihedrals]
            The dihedral angles, in radians
        """
        # extract the dihedrals from trajectory X
        return dihedrals


##############################################################################
# Two different patterns for "meta" vectorizers.

# One is a "merge" that operates one or more vectorizers in parallel and then
# stacks the results together, to get a single set of features that includes
# the features calculated by each of the vectorizers.

# The second is a "pipeline" that sends the data through the vectorizers in
# a sequence, like first calculating distances, using those as the input for,
# say, PCA, and then using the transformed output from PCA as input to tICA.
##############################################################################


class MergingTransformer(BaseEstimator, TransformerMixin):
    """Transforms a trajectory by applying a collection of other vectorizers,
    stacking the results together.

    This transformer applies a series of other transformers to the input
    trajectories, giving a single set of features for each frame which might
    each be computed by a different vectorizer. For example, a
    MergingVectorizer built from a DistanceVectorizer and AngleFeaturizer
    would, when `transform`ing a trajectory, return both the calculated
    distances and angles in a single array.

    Examples
    --------
    >>> a = AngleVectorizer([[0, 1, 2]])
    >>> b = DistanceVectorizer([[1,2]])
    >>> X = md.load('trajectory.h5')
    >>> features = MergingTransformer([a, b]).fit_transform(X)
    """

    def __init__(self, transformers):
        self.transformers = transformers

    def fit(self, X):
        """Fit this vectorizer to training data

        Returns
        -------
        self
        """
        for transformer in self.transformers:
            transformer.fit(X)
        return self

    def transform(self, X):
        """Extract features from each frame in a trajectory

        X_new : numpy array of shape [n_frames, features]
            The features for each frame, computed by applying each of the
            vectorizers and concatenating their results
        """
        return np.hstack([v.transform(X) for v in self.transformers])


class PipelineTransformer(BaseEstimator, TransformerMixin):
    """Transforms a trajectory by applying a pipeline sequence of
    transformations in order, with the results of one feeding the input
    to the subsequent transformer.
    
    Examples
    --------
    >>> a = DistanceVectorizer([[1,2]])
    >>> b = tICA(n_components=2)

    >>> X = md.load('trajectory.h5')
    >>> features = PipelineTransformer([a, b]).fit_transform(X)
    """
    
    def __init__(self, transformers):
        self.transformers = transformers
    
    def fit(self, X):
        """Fit this vectorizer to training data

        Returns
        -------
        self
        """
        for transformer in self.transformers:
            X = transformer.fit_transform(X)
        return self
    
    def transform(self, X):
        """Apply this sequence of transformations to new data
        
        Returns
        -------
        X_new : numpy array of shape [n_frames, features]
            The features for each frame, computed by the chain of transformers
            operating sequentially
        """

        for transformer in self.transformers:
            X = transformer.transform(X)
        return X
    
    def fit_transform(self, X):
        """Fit and apply this sequence of transformations

        Returns
        -------
        X_new : numpy array of shape [n_frames, features]
            The features for each frame, computed by the chain of transformers
            operating sequentially
        """
        for transformer in self.transformers:
            X = transformer.fit_transform(X)
        return X
        

##############################################################################
# tICA is also a vectorizer
##############################################################################


class tICA(BaseEstimator, TransformerMixin):
    """Time-structure based independent component analysis (tICA)

    Linear dimensionality reduction of the data, keeping the most slowly
    decorrelating components of the data to project the data into a lower
    dimensional space


    Parameters
    ----------
    n_components : int, None
        Number of components to keep.
        if n_components is not set all components are kept::

            n_components == min(n_samples, n_features)

    """

    def __init_(self, n_components):
        self.n_components = n_components

    def fit(self, X):
        """Fit the model with one or more vectorized trajectories.
        """
        # ... train tICA
        self.mean_ = do the means need to be subtracted out?
        self.components_ = the retained tICA eigenvectors

    def transform(self, X):
        """Apply the dimensionality reduction on a vectorized trajectory

        Parameters
        ----------
        X : numpy array of shape [n_frames, n_features]
            New data, where n_frames is the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : anumpy array of shape [n_frames, n_components]
        """
        # maybe center X ?
        X_transformed = np.dot(X, self.components_.T)
        return X_transformed
