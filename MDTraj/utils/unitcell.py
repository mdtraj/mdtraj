##############################################################################
# Imports
##############################################################################

import warnings
import numpy as np

##############################################################################
# Functions
##############################################################################


def lengths_and_angles_to_box_vectors(a_length, b_length, c_length, alpha, beta, gamma):
    """Convert from the lengths/angles of the unit cell to the box
    vectors (Bravais vectors). The angles should be in degrees.

    Parameters
    ----------
    a_length : scalar or np.ndarray
        length of Bravais unit vector **a**
    b_length : scalar or np.ndarray
        length of Bravais unit vector **b**
    c_length : scalar or np.ndarray
        length of Bravais unit vector **c**
    alpha : scalar or np.ndarray
        angle between vectors **b** and **c**, in degrees.
    beta : scalar or np.ndarray
        angle between vectors **c** and **a**, in degrees.
    gamma : scalar or np.ndarray
        angle between vectors **a** and **b**, in degrees.

    Returns
    -------
    a : np.ndarray
        If the inputs are scalar, the vectors will one dimesninoal (length 3).
        If the inputs are one dimension, shape=(n_frames, ), then the output
        will be (n_frames, 3)
    b : np.ndarray
        If the inputs are scalar, the vectors will one dimesninoal (length 3).
        If the inputs are one dimension, shape=(n_frames, ), then the output
        will be (n_frames, 3)
    c : np.ndarray
        If the inputs are scalar, the vectors will one dimesninoal (length 3).
        If the inputs are one dimension, shape=(n_frames, ), then the output
        will be (n_frames, 3)

    Examples
    --------
    >>> lengths_and_angles_to_box_vectors(1, 1, 1, 90.0, 90.0, 90.0)
    (array([1, 0, 0]), array([  6.12323400e-17,   1.00000000e+00,   0.00000000e+00]), array([  6.12323400e-17,   6.12323400e-17,   1.00000000e+00]))


    Notes
    -----
    This code is adapted from gyroid, which is licensed under the BSD
    http://pythonhosted.org/gyroid/_modules/gyroid/unitcell.html
    """
    if np.all(alpha < 2*np.pi) and np.all(beta < 2*np.pi) and np.all(gamma < 2*np.pi):
        warnings.warn('All your angles were less than 2*pi. Did you accidently give me radians?')

    alpha = alpha * np.pi / 180
    beta = beta * np.pi / 180
    gamma = gamma * np.pi / 180

    a = np.array([a_length, np.zeros_like(a_length), np.zeros_like(a_length)])
    b = np.array([b_length*np.cos(gamma), b_length*np.sin(gamma), np.zeros_like(b_length)])
    cx = c_length*np.cos(beta)
    cy = c_length*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))
    cz = np.sqrt(c_length*c_length - cx*cx - cy*cy)
    c = np.array([cx,cy,cz])

    if not a.shape == b.shape == c.shape:
        raise TypeError('Shape is messed up.')

    return a.T, b.T, c.T


def box_vectors_to_lengths_and_angles(a, b, c):
    """Convert box vectors into the lengths and angles definining the box

    Parameters
    ----------
    a : np.ndarray
        the vector defining the first edge of the periodic box (length 3), or
        an array of this vector in multiple frames, where a[i,:] gives the
        length 3 array of vector a in each frame of a simulation
    b : np.ndarray
        the vector defining the second edge of the periodic box (length 3), or
        an array of this vector in multiple frames, where b[i,:] gives the
        length 3 array of vector a in each frame of a simulation
    c : np.ndarray
        the vector defining the third edge of the periodic box (length 3), or
        an array of this vector in multiple frames, where c[i,:] gives the
        length 3 array of vector a in each frame of a simulation

    Examples
    --------
    >>> a = np.array([2,0,0], dtype=float)
    >>> b = np.array([0,1,0], dtype=float)
    >>> c = np.array([0,1,1], dtype=float)
    >>> l1, l2, l3, alpha, beta, gamma = box_vectors_to_lengths_and_angles(a, b, c)
    >>> (l1 == 2.0) and (l2 == 1.0) and (l3 == np.sqrt(2))
    True
    >>> np.abs(alpha - 45) < 1e-6
    True
    >>> np.abs(beta - 90.0) < 1e-6
    True
    >>> np.abs(gamma - 90.0) < 1e-6
    True

    Returns
    ----------
    a_length : scalar or np.ndarray
        length of Bravais unit vector **a**
    b_length : scalar or np.ndarray
        length of Bravais unit vector **b**
    c_length : scalar or np.ndarray
        length of Bravais unit vector **c**
    alpha : scalar or np.ndarray
        angle between vectors **b** and **c**, in degrees.
    beta : scalar or np.ndarray
        angle between vectors **c** and **a**, in degrees.
    gamma : scalar or np.ndarray
        angle between vectors **a** and **b**, in degrees.
    """
    if not a.shape == b.shape == c.shape:
        raise TypeError('Shape is messed up.')
    if not a.shape[-1] == 3:
        raise TypeError('The last dimension must be length 3')
    if not (a.ndim in [1,2]):
        raise ValueError('vectors must be 1d or 2d (for a vectorized '
                         'operation on multiple frames)')
    last_dim = a.ndim-1

    a_length = np.sqrt(np.sum(a*a, axis=last_dim))
    b_length = np.sqrt(np.sum(b*b, axis=last_dim))
    c_length = np.sqrt(np.sum(c*c, axis=last_dim))

    # we allow 2d input, where the first dimension is the frame index
    # so we want to do the dot product only over the last dimension
    alpha = np.arccos(np.einsum('...i, ...i', b, c) / (b_length * c_length))
    beta = np.arccos(np.einsum('...i, ...i', c, a) / (c_length * a_length))
    gamma = np.arccos(np.einsum('...i, ...i', a, b) / (a_length * b_length))

    # convert to degrees
    alpha = alpha * 180.0 / np.pi
    beta = beta * 180.0 / np.pi
    gamma = gamma * 180.0 / np.pi

    return a_length, b_length, c_length, alpha, beta, gamma
