from __future__ import print_function, division
import sysconfig
from distutils.ccompiler import new_compiler
import cffi
import numpy as np

from .six import iteritems, PY3

__all__ = ['cdata', 'find_library']

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

    def __init__(self):
        ffi = cffi.FFI()
        # some platforms dont have all of these, especially wierd compilers or 32 bit machines
        nptypes = [getattr(np, name) for name in self.nptype_names if hasattr(np, name)]
        nptype_descr = {'%s%d' % (dtype.kind, dtype.itemsize): dtype for dtype in map(np.dtype, nptypes)}

        casts = {}
        for code, names in iteritems(self.ctypes):
            for name in names: 
                casts[nptype_descr['%s%d' % (code, ffi.sizeof(name))]] = name + ' *'
        # casts is a dict that helps us cast numpy arrays, like
        # {np.float32 : 'float *', np.int32: 'int *'}
        
        
        self._ffi = ffi
        self._casts = casts

    def __call__(self, ndarray):
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
