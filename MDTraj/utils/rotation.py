##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: John D. Chodera and Kyle A. Beauchamp, Kim Branson,
#          Imran Haque, Michael Shirts
# Contributors: Robert McGibbon
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
#
# This code is adapted from Yank: GPU-accelerated calculation of ligand
# binding affinities in implicit and explicit solvent using alchemical
# free energy methodologies, copyright 2009-2011 University of California,
# Berkeley, Vertex Pharmaceuticals, Stanford University, University of Virginia,
# and the Authors. Yank is licensed under the GPU Lesser General Public
# License.
##############################################################################

import numpy as np
from mdtraj.utils import check_random_state

__all__ = ["rotation_matrix_from_quaternion", "uniform_quaternion"]



def rotation_matrix_from_quaternion(q):
    """Compute a 3x3 rotation matrix from a given quaternion (4-vector).
    
    Parameters
    ----------
    q : np.ndarray, size=(..., 4)
        Quaternion or array of quaternions (need not be normalized, zero norm OK)

    Returns
    -------
    Rq : np.ndarray, size=(..., 3, 3)
        Orthogonal rotation matrices corresponding to quaternion(s) q. Given,
        for example, an input shape of (m, 4), Rq will have shape (m, 3, 3).

    Examples
    --------
    >>> q = np.array([0.1, 0.2, 0.3, -0.4])
    >>> Rq = rotation_matrix_from_quaternion(q)

    See Also
    --------
    uniform_quaternion

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
    """

    w = q[..., 0]; x = q[..., 1]; y = q[..., 2]; z = q[..., 3]
    # w, x, y, z = q
    Nq = np.sum(q**2, axis=q.ndim-1)

    if q.shape == (4,):
        s = 2.0 / Nq if Nq > 0 else 0.0
    else:
        with np.errstate(divide='ignore'):
            s = 2.0 / Nq
            s[Nq <= 0.0] = 0.0

    X = x*s; Y = y*s; Z = z*s
    wX = w*X; wY = w*Y; wZ = w*Z
    xX = x*X; xY = x*Y; xZ = x*Z
    yY = y*Y; yZ = y*Z; zZ = z*Z
    Rq = np.array([[ 1.0-(yY+zZ),       xY-wZ,        xZ+wY  ],
                   [      xY+wZ,   1.0-(xX+zZ),       yZ-wX  ],
                   [      xZ-wY,        yZ+wX,   1.0-(xX+yY) ]])

    if q.shape == (4,):
        return Rq
    return np.rollaxis(Rq, Rq.ndim-1)


def uniform_quaternion(size=None, random_state=None):
    """Generate uniform normalized quaternion 4-vectors

    Parameters
    ----------
    size : int or tuple of ints, optional
        Defines the shape of the returned array of quaternions. If None
        (the default), returns a quaternion 4-vector.
   random_state : integer or numpy.RandomState, optional
        The generator used for random numbers. If an integer is given,
        it fixes the seed. Defaults to the global numpy random number
        generator.

    Returns
    -------
    out : ndarray
        Array of quaternion 4-vectors. Given a ``size`` of, for example,
        ``(m,n)``, ``m*n`` samples are generated.  Because each sample is
        `4-dimensional, the output shape is ``(m,n,4)``. If no shape is
        specified, a single (`4`-D) sample is returned.

    See Also
    --------
    rotation_matrix_from_quaternion

    References
    ----------
    .. [1] K. Shoemake. Uniform random rotations. In D. Kirk, editor, Graphics
       Gems III, pages 124-132. Academic, New York, 1992.
    .. [2] http://planning.cs.uiuc.edu/node198.html
    """
    random = check_random_state(random_state)
    if size is None:
       randsize = (3, 1)
    elif isinstance(size, tuple):
        size = tuple(int(i) for i in size)
        randsize = (3,) + size
    else:
        randsize = (3, int(size))

    u = random.random_sample(randsize)
    q = np.array([np.sqrt(1 - u[0]) * np.sin(2*np.pi * u[1]),
                  np.sqrt(1 - u[0]) * np.cos(2*np.pi * u[1]),
                  np.sqrt(u[0])     * np.sin(2*np.pi * u[2]),
                  np.sqrt(u[0])     * np.cos(2*np.pi * u[2])])

    if size is None:
        return q[:, 0]
    return np.rollaxis(q, 0, q.ndim)
