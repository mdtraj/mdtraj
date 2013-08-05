from distutils.ccompiler import new_compiler
import cffi
import numpy as np

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
    nptypes = [np.float16, np.float32, np.float, np.float64, np.double, np.float128,
               np.int16, np.int32, np.int, np.int64, np.long,
               np.uint16, np.uint32, np.uint64, np.uint]

    def __init__(self):
        ffi = cffi.FFI()
        nptype_descr = {'%s%d' % (dtype.kind, dtype.itemsize): dtype for dtype in map(np.dtype, self.nptypes)}

        casts = {}
        for code, names in self.ctypes.iteritems():
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
    if not hasattr(path, '__iter__'):
        path = [path]
    return compiler.find_library_file(path, name)