##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors:
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
##############################################################################

from __future__ import print_function, division
from distutils.ccompiler import new_compiler
try:
    import sysconfig  # py3
except ImportError:
    from distutils import sysconfig

import numpy as np
from mdtraj.utils.six import iteritems, PY3

__all__ = ['cpointer', 'find_library']

class cpointer(object):
    """Function that extracts a cffi pointer from a numpy array. The function
    is actually structured as an instance of a callable class so that some
    initialization -- building a table of the appropriate numpy dtype ->
    native c casts can be constucted at initialization.
    """
    ctypes = {
        'f': ['float', 'double', 'long double'],
        'i': ['short', 'int', 'long'],
        'u': ['unsigned short', 'unsigned int', 'unsigned long']
        }
    nptype_names = ['float16', 'float32', 'float', 'float64', 'double', 'float128',
                   'int16', 'int32', 'int', 'int64', 'long',
                   'uint16', 'uint32', 'uint64', 'uint']

    def _delayed_init(self):
        import cffi
        ffi = cffi.FFI()
        # some platforms dont have all of these, especially wierd compilers or 32 bit machines
        nptypes = [getattr(np, name) for name in self.nptype_names if hasattr(np, name)]
        nptype_descr = dict([('%s%d' % (dtype.kind, dtype.itemsize), dtype) for dtype in map(np.dtype, nptypes)])

        casts = {}
        for code, names in iteritems(self.ctypes):
            for name in names:
                casts[nptype_descr['%s%d' % (code, ffi.sizeof(name))]] = name + ' *'
        # casts is a dict that helps us cast numpy arrays, like
        # {np.float32 : 'float *', np.int32: 'int *'}
        self._ffi = ffi
        self._casts = casts

    def __call__(self, ndarray):
        if not hasattr(self, '_ffi'):
            self._delayed_init()
        if not isinstance(ndarray, np.ndarray):
            raise TypeError('ndarray must be an instance of numpy.ndarray')
        if not ndarray.flags.contiguous:
            raise ValueError('ndarray must be contiguous')
        return self._ffi.cast(self._casts[ndarray.dtype], ndarray.ctypes.data)

# overshadow the class with an instance of the class -- which is just effectively
# a callable function w/ some precomputed data
cpointer = cpointer()


def find_library(path, name):
    """Find a shared library
    """
    compiler = new_compiler()
    if isinstance(path, str):
        path = [path]
    names = [name]

    soabi = sysconfig.get_config_var('SOABI')
    if soabi is not None:
        names.append('%s.%s' % (names[0], soabi))

    for name in names:
        result = compiler.find_library_file(path, name)
        if result is not None:
            return result

    return None
