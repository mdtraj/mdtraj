
from mdtraj.vectorizer import (BaseModeller, EstimatorMixin,
                               TransformerMixin, UpdateableEstimatorMixin)
import numpy as np

class tICAVectorizer(BaseModeller, UpdateableEstimatorMixin, TransformerMixin):
    """
    Class for doing time-structure based Independent Component Analysis
    (tICA). See (Schwantes, CR and Pande, VS. JCTC, 2013, 9 (4), pp 2000-09)
    for further details.

    Briefly, tICA finds the linear combinations that maximize the projected
    autocorrelation function. This means that the tICA components correspond
    to the slowest decorrelating linear combinations of the input coordinates.
    """
    def __init__(self, lag, n_components=None, pca_cutoff=1E-8):
        """
        initialize the tICAVectorizer object

        Parameters
        ----------
        lag : int
            lag time to use in calcualting the timelag correlation matrix
            the unit is in frames
        n_components : int, optional
            number of components to project onto. You can set this later...
        pca_cutoff : float
            The estimated covariance matrix must be strictly positive definite.
            To ensure this, we can work in the PCA space defined by all PCs
            with nonzero variance. The pca_cutoff is the cutoff for defining
            zero variance (Default: 1E-8)
        """
        
        self.n_components = n_components
        self.lag = lag

        # set up containers for running sums
        self.running_corr_0_dt_ = None
        self.running_corr_0_0_ = None
        self.running_corr_dt_dt_ = None
        self.running_sum_0_ = None
        self.running_sum_dt_ = None

        self.total_frames_ = 0

        self._have_estimate_ = False

        
    def fit_update(self, X):
        """
        Update the internal state with new data, X

        Parameters
        ----------
        X : np.ndarray or list of np.ndarrays
            Data to add to the estimate of the tICA matrices. This can be a list
            of numpy arrays, or a single numpy array. Each array should be two-
            dimensional: (n_samples, n_coordinates). 

        Returns
        -------
        self
        """

        if not isinstance(X, list):
            X = [X]

        self._have_estimate_ = False
        # we have updated the data, so we no longer have the PCs.

        for row in X:

            if not isinstance(row, np.ndarray):
                raise RuntimeError("data must be numpy.ndarray's or a list of arrays")

            shape = row.shape

            if len(shape) == 1 or shape[0] <= self.lag:
                logger.warn("row is too short, not using this data")

            if len(shape) > 2:
                raise RuntimeError("data cannot be more than two-dimensional")

            n_features = row.shape[1]

            if not self.cov_mat_ is None:
                self.running_corr_0_0_ = np.zeros((n_features, n_features))
                self.running_corr_0_dt_ = np.zeros((n_features, n_features))
                self.running_corr_dt_dt_ = np.zeros((n_features, n_features))
                self.running_sum_0_ = np.zeros(n_features)
                self.running_sum_dt_ = np.zeros(n_features)

            elif n_features != self.cov_mat_.shape[0]:
                raise RuntimeError("data does not match the shape of the internal state.")

            self.running_corr_mat_0_0_ += row[:-self.lag].T.dot(row[:-self.lag])
            self.running_corr_mat_dt_dt_ += row[self.lag:].T.dot(row[self.lag:])
            self.running_corr_mat_0_dt_ += row[:-self.lag].T.dot(row[self.lag:])
            self.running_sum_0_ += row[:-self.lag].sum(0)
            self.running_sum_dt_ += row[self.lag:].sum(0)

            self.total_samples_ += row.shape[0] - self.lag

        return self


    def fit(self, X):
        """
        calculate the slowest components for data, X

        Parameters
        ----------
        X : np.ndarray or list of np.ndarrays
            data or list of numpy arrays
        
        Returns
        -------
        self
        """
        self.clear()

        if isinstance(X, list):
            for row in X:
                self.fit_update(row)

        else:
            self.fit_update(row)

        self.compute_components()

        return self

        
    def clear(self):
        """
        clear the internal state, to analyze new data with tICA
        """
        
        super(self, tICAVectorizer).clear()

        self.total_frames_ = 0
        self._have_estimate_ = False


    def compute_components(self):
        """
        compute the solutions to the tICA problem:

        C v = w S v

        We will restrict our solutions to linear combinations of
        non-zero variance principal components. As long as we really are 
        ignoring only the zero-variance PCs this will give the right
        solutions.

        See discussion in Shukla, D et. al. In preperation. Or email
        Christian Schwantes (schwancr@stanford.edu) for a detailed
        explanation.
        """
        
        # first we have to do PCA, since odds are our covariance matrix
        # is not positive definite

        self.mean_ = (self.running_sum_0_ + self.running_sum_dt_) / (2. * float(self.total_frames_))
        outer_mean = np.outer(self.mean_, self.mean_)
        cov_mat = (self.running_corr_0_0_ + self.running_corr_dt_dt_) / (2. * float(self.total_frames_))
        cov_mat = cov_mat - outer_mean

        timelag_corr_mat = (self.running_corr_0_dt_ + self.running_corr_0_dt_.T) / (2. * float(self.total_frames_))
        timelag_corr_mat = timelag_corr_mat - outer_mean

        pca_vals, pca_vecs = np.linalg.eigh(cov_mat)

        ind = np.where(pca_vals > self.pca_cutoff)[0]
        pca_vecs = pca_vecs[:, ind]

        lhs = pca_vecs.T.dot(timelag_corr_mat).dot(pca_vecs)
        rhs = pca_vecs.T.dot(cov_mat).dot(pca_vecs)

        vals, vecs = scipy.linalg.eig(lhs, b=rhs)

        if vals.imag.max() > 1E-10:
            logger.warn("there are non-real eigenvalues in this solution. "
                        "You should probably use a larger pca_cutoff.")

        else:
            vals = vals.real
            vecs = vecs.real

        ind = np.argsort(vals)[::-1]

        self.vals_ = vals[ind]
        self.vecs_ = vecs[:, ind]

        self._have_estimate_ = True


    def transform(self, X):
        """
        Transform some data, X, onto the slowest n_components

        Parameters
        ----------
        X : np.ndarray or list of np.ndarray's
            data to project onto the top n_components. Should be a single two-
            dimensional array (n_samples, n_coordinates) or a list of arrays
            
        Returns
        -------
        proj_X : np.ndarray or list of np.ndarray's
            projected data

        """

        if self.n_components is None:
            raise RuntimeError("need to set n_components")

        if not self._have_estimate_:
            self.compute_components()

        if isinstance(X, list):
            return_list = True

        else:
            X = [X]
            return_list = False

        top_tics = self.vecs_[:, :self.n_components]

        proj_X = []
        for row in X:
            if not isinstance(row, np.ndarray):
                raise RuntimeError("data contains rows that are not np.ndarray's")

            shape = row.shape
            if len(shape) == 1:
                row = row.reshape((1, -1))

            if len(shape) > 2:
                raise RuntimeError("data cannot be more than two-dimensional")

            n_features = row.shape[1]

            if n_features != top_pcs.shape[0]:
                raise RuntimeError("data is not the right shape")

            proj_X.append(row.dot(top_tics))
            # are you supposed to subtract the mean before projecting?
            # if so, then this is the correct line:

            # proj_X.append((row - self.mean_).dot(top_pcs)

            # but, this just adds a constant vector to each point, (-self.mean_.dot(top_pcs))
            # so I don't think it actually matters..

        if return_list:
            return proj_X
        else:
            return proj_X[0]

