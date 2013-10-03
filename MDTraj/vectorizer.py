import inspect
import itertools
import tables
import numpy as np
import mdtraj as md


class BaseModeller(object):
    """Base class for all statistical modelling

    Notes
    -----
    All subclasses should specified all their parameters
    at the cass level in their __init__ as explicit keyword
    args (no *args, **kwargs)
    """

    def get_params(self):
        """Get parameters for this modeller

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr(self, key)
            out[key] = value

        return out

    def set_params(self, **params):
        """Set the parameters of this modeller

        Returns
        -------
        self
        """
        if not params:
            return self

        valid_params = self._get_param_names()
        for key, value in params.iteritems():
            if key not in valid_params:
                raise ValueError('Invalid parameter %s for modeller %s'
                                 % (key, self.__class__.__name__))
            setattr(self, key, value)

        return self

    @classmethod
    def _get_param_names(cls):
        """Get the parameter names for this modeller"""
        try:
            init = cls.__init__
            args, varargs, kw, default = inspect.getargspec(init)
            if varargs is not None:
                raise RuntimeError("mdtraj modellers should always "
                                   "specify their parameters in the signature"
                                   " of their __init__ (no varargs)."
                                   " %s doesn't follow this convention."
                                   % (cls, ))
            # remove self
            args.pop(0)
        except TypeError:
            # no explicit __init__
            args = []
        
        args.sort()
        return args


class TransformerMixin(object):
    """"Mixin class for all transformers"""

    def transform(self, X):
        """Transform a dataset X from one represenation/basis to another

        Parameters
        -----------
        X : numpy array, trajectory, list of numpy arrays, or list of trajectories
            If X is an individual array or trajectory, we transform it
            into the new space indivually. If X is a list of arrays or
            trajectories, the whole dataset is transformed, and we return a
            list of new arrays

        Returns
        -------
        X_new : numpy array of shape [n_samples, n_features_new], or list of such arrays
            Transformed data, a represenation of the input data X in the new
            space produced by this transformer
        """
        raise NotImplementedError()


class EstimatorMixin(object):
    """Mixin class for all estimators"""

    def fit(self, X):
        """Fit this vectorizer to training data

        Parameters
        ----------
        X : any
            The dataset to fit the estimator with. This dataset should
            encompas all of the data. Repeated calls to `fit` do not
            incrementally update the estimator.

        Returns
        -------
        self
        """
        return self

    def _get_estimate_names(self):
        """Each estimator object, when `fit` on a dataset, produces
        estimates which are instance variables on the class whose names
        end with an underscore, e.g. self.mean_ or self.variance_.

        Returns
        -------
        names : list of strings
             A list of the names of each of the variables currently
             estimated by this estimator.
        """
        # I'm not sure how robust this is...
        return [e for e in self.__dict__.keys() if e.endswith('_') and not e.startswith('_')]

    def to_pytables(self, parentnode):
        """Serialize this estimator to a PyTables group, attaching
        it to the parent node.

        It this estimator wraps a series of other estimators, as
        in a pipeline, it should call to_pytables on its children
        recursively.

        Parameters
        ----------
        parentnode : tables.Group
            The parent node in the HDF5 hierachy that this group
            should be attached to.
        """
        group = tables.Group(parentnode, name=self.__class__.__name__, new=True)

        # Supported types are int, float, str and numpy arrays. Atomic types (int, float, str)
        # will go in a Table named params, and numpy arrays will each go in an individual
        # Array

        params = self.get_params()
        estimates = {k: getattr(self, k) for k in self._get_estimate_names()}

        simple_typemap = {int: tables.Int64Col(), float: tables.FloatCol(), str: tables.StringCol(1024),
                          np.int64: tables.Int64Col(), np.int32: tables.Int32Col(), np.float32: tables.Float32Col(),
                          np.float64: tables.Float64Col(), bool: tables.BoolCol}

        table_description = {}
        table_entries = {}
        for key, value in itertools.chain(params.iteritems(), estimates.iteritems()):
            if value is None:
                # just skip Nones. They're indicated by their absense.
                continue
            if isinstance(value, np.ndarray):
                parentnode._v_file.create_array(group, key, obj=value)
            else:
                try:
                    table_description[key] = simple_typemap[type(value)]
                    table_entries[key] = value
                except KeyError:
                    raise RuntimeError("I don't know how to serialize the parameter %s (type=%s) to pytables" 
                                       % (key, type(value)))

        table = parentnode._v_file.create_table(group, 'params', expectedrows=1, description=table_description)
        for key, value in table_entries.iteritems():
            table.row[key] = value
        
        return group

    @classmethod
    def from_pytables(cls, group):
        """Instantiate a copy of this estimator from a pytables group.
        This performs the inverse of to_pytables.
        """

        instance = cls()
        # instance.set_params(extract parameters from group)
        # if this group contains any nested subgroups and cls
        # is some kind of meta-modeller like a pipeline, then
        # this method needs to recursively call from_pytables
        # on the subgroups, to fully instantiate instance.
        return instance


if __name__ == '__main__':
    class Test(BaseModeller, EstimatorMixin):
        def __init__(self, a=1, b='a', c=None):
            self.a = 1
            self.b = b
            self.c = c
            
        def fit(self, X):
            self.mean_ = np.array([1,2,3])
            self.mean0_ = self.mean_[1] 
            self.mean1_ = 4.0
            return self

    t = Test().fit(None)
    f = tables.open_file('f.h5', 'w')
    t.to_pytables(f.root)
    print f
    print f.root.Test.params.row[:]

    f.close()



class UpdatableEstimatorMixin(EstimatorMixin):
    def fit_update(self, X):
        """Update the statistical model described by this estimator by
        exposing it to new data, without "forgetting" the data that
        it has seen in any previous calls to `fit` or `fit_update`

        Parameters
        ----------
        X : any
            The dataset to update the estimator with.

        Returns
        -------
        self
        """
        raise NotImplementedError()

    def fit(self, X):
        """Fit this estimator to training data

        Parameters
        ----------
        X : any
            The dataset to fit the estimator with. This dataset should
            encompas all of the data. Repeated calls to `fit` do not
            incrementally update the estimator.

        Returns
        -------
        self
        """
        self.clear()
        if isinstance(X, list):
            for x in X:
                self.fit_update(x)
        else:
            self.fit_update(X)
        return self

    def clear(self):
        """Clear the state of this estimator, so that it can be refit
        on fresh data without retaining knowledge of data it has been
        previously exposed to.

        All instance variables ending in '_' (estimated quantities)
        will be set to None
        """

        for name in self.__dict__.keys():
            if name.endswith('_'):
                setattr(self, name, None)


##############################################################################
# Vectorizers: Take a trajecory as input, and return a [n_frames, n_features]
# array
##############################################################################


class PositionVectorizer(BaseModeller, TransformerMixin):
    """Transforms a molecular dynamics trajectory into a multvariate
    timeseries of the positions of specified atoms in cartesian space

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
        self.alignment_indices = alignment_indices
        self.reference = reference

    @property
    def _target(self):
        """_target is the cartesian coordinates of the alignment
        atoms of the first frame of the reference trajectory."""
        if not hasattr(self, '__target') or not self.__target:
            self.__target = self._reference.xyz[0, self.alignment_indices]
        return self.__target

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
        if isinstance(X, list):
            return map(self._transform, X)
        return self._transform(X)

    def _transform(self, X):
        # This can be probably done more efficiently in C
        X_new = np.empty_like(X.xyz)
        for i in range(X.n_frames):
            X_new[i] = md.geometry.alignment.transform(X.xyz[i, self.alignment_indices], self._target)

        return X_new


class DistanceVectorizer(BaseModeller, TransformerMixin):
    """Transforms a molecular dynamics trajectory into a multvariate 
    timeseris of pairwise distances between specified atoms

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
    >>> distances = DistanceVectorizer([[1,2]]).transform(X)
    """

    def __init__(self, pair_indices, use_periodic_boundries=False):
        self.pair_indices = pair_indices
        self.use_periodic_boundries = use_periodic_boundries

    def transform(self, X):
        """Extract the distance between pairs of atoms from a trajectory

        Parameters
        ----------
        X : Trajectory, or list of Trajectories
            One or more molecular dynamics trajectories

        Returns
        -------
        d : numpy array of shape [n_frames, n_distances]
            One or more arrays of pairwise distances
        """
        if isinstance(X, list):
            return map(self._transform, X)
        return self._transform(X)

    def _transform(self, X):
        return md.geometry.compute_distances(X, self.pair_indices, self.use_periodic_boundaries)


class AngleVectorizer(BaseModeller, TransformerMixin):
    """Transforms a molecular dynamics trajectory into a multivariate
    timeseries of the angles between specific atoms

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
    >>> distances = AngleVectorizer([[0, 1, 2]]).transform(X)
    """

    def __init__(self, triplet_indices):
        self.triplet_indices = triplet_indices

    def transform(self, X):
        """Extract the angles between atoms from a trajectory

        Parameters
        ----------
        X : Trajectory or list of Trajectories
            One or more molecular dynamics trajectories

        Returns
        -------
        d : numpy array of shape [n_frames, n_angles]
            One or more arrays of angles, in radians
        """
        if isinstance(X, list):
            return map(self._transform, X)
        return self._transform(X)
    
    def _transform(self, X):
        return md.geometry.compute_angles(X, self.triplet_indices)


class DihedralVectorizer(BaseModeller, TransformerMixin):
    """Transforms a trajectory into a multivariate timeseries of the
    torsion angles between specific atoms

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
    >>> distances = DihedralVectorizer([[0, 1, 2, 3]]).transform(X)
    """
    def __init__(self, quartet_indices):
        self.quartet_indices = quartet_indices

    def transform(self, X):
        """Extract torsion angles from a trajectory

        Parameters
        ----------
        X : Trajectory or list of trajectories
            One or more molecular dynamics trajectories

        Returns
        -------
        d : numpy array of shape [n_frames, n_dihedrals]
            One or more arrays of dihedral angles, in radians
        """
        if isinstance(X, list):
            return map(self._transform, X)
        return self._transform(X)

    def _transform(self, X):
        return md.geometry.compute_dihedrals(X, self.quartet_indices)


##############################################################################
# Two different patterns for "meta" vectorizers.

# One is a "merge" that operates one or more vectorizers in parallel and then
# stacks the results together, to get a single set of features that includes
# the features calculated by each of the vectorizers.

# The second is a "pipeline" that sends the data through the vectorizers in
# a sequence, like first calculating distances, using those as the input for,
# say, PCA, and then using the transformed output from PCA as input to tICA.
##############################################################################


class MergingTransformer(BaseModeller, TransformerMixin):
    """Transforms a trajectory or timeseries by applying a collection of
    other transformers and stacking the results together.

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
    >>> features = MergingTransformer([a, b]).transform(X)
    """

    def __init__(self, transformers):
        self.transformers = transformers

    def transform(self, X):
        """Extract features from each frame in a trajectory

        X_new : numpy array of shape [n_frames, features]
            The features for each frame, computed by applying each of the
            vectorizers and concatenating their results
        """
        return np.hstack([v.transform(X) for v in self.transformers])


class PipelineTransformer(BaseModeller, TransformerMixin):
    """Transforms a trajectory or timeseries by applying a pipeline
    sequence of transformations in order, with the results of one
    feeding the input to the subsequent transformer.

    Examples
    --------
    >>> a = DistanceVectorizer([[1,2]])
    >>> b = tICA(n_components=2)

    >>> X = md.load('trajectory.h5')
    >>> features = PipelineTransformer([a, b]).transform(X)
    """

    def __init__(self, transformers):
        self.transformers = transformers

    def transform(self, X):
        """Apply this sequence of transformations to new data

        Parameters
        ----------
        X : any
            One or more trajectories or arrays, suitable as input for the
            first transformer in the pipeline

        Returns
        -------
        X_new : numpy array of shape [n_frames, features]
            The features for each frame, computed by the chain of transformers
            operating sequentially
        """

        for t in self.transformers:
            X = t.transform(X)
        return X


##############################################################################
# tICA  is both a transformer and an estimator
##############################################################################


class tICA(BaseModeller, TransformerMixin, UpdatableEstimatorMixin):
    """Time-structure based independent component analysis (tICA)

    Linear dimensionality reduction of multivariate timeseries data, keeping
    the most slowly decorrelating components of the data, with which the
    timeseries can be projected into a lower dimensional space

    Parameters
    ----------
    n_components : int, None
        Number of components to keep.
        if n_components is not set all components are kept::

            n_components == min(n_samples, n_features)

    """

    def __init_(self, n_components):
        self.n_components = n_components

        # estimated quantites
        self.components_ = None
        self.eigenvalues_ = None
        self.means_ = None

        self.covariance_ = None
        self.time_lag_correlation_ = None

    def fit_update(self, X):
        """Update the tICA estimator with one or more featurized trajectories

        Parameters
        ----------
        X : one or more numpy array of shape [n_frames, n_features]
             One (or more) multivariate timeseries from the input dataset

        Returns
        -------
        self
        """
        # update self.covariance_ and self.time_lag_correlation_
        pass

    def compute_components(self):
        """Compute the slow components of the dataset

        This operation requires diagonalizing a matrix of size
        [n_features, n_features], so it can be somewhat costly
        """
        pass

    def transform(self, X):
        """Apply the tICA dimensionality reduction to a multivariate timeseries

        Parameters
        ----------
        X : numpy array of shape [n_frames, n_features]
            New data, where n_frames is the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : anumpy array of shape [n_frames, n_components]
        """
        if self.covariance_ is None or self.time_lag_correlation_ is None:
            self.compute_components()

        # maybe center X ?
        X_transformed = np.dot(X, self.components_.T)
        return X_transformed
