__all__ = ["rg", "angle", "distance", "dihedral", "internal", "alignment"]

import os
import warnings
try:
    import cffi
    from mdtraj.utils.ffi import find_library
    _HAVE_OPT = None   # not sure if we have the library yet
except ImportError:
    warnings.warn('Optimized distance library requires the "cffi" package, '
                  'which is installable with easy_install or pip via '
                  '"pip install cffi" or "easy_install cffi".')
    _HAVE_OPT = False  # we definitely don't have the library

if _HAVE_OPT is not False:
    # lets try to open the library
    ffi = cffi.FFI()
    ffi.cdef('''int dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
                             float* distance_out, float* displacement_out,
                             const int n_frames, const int n_atoms, const int n_pairs);''')
    ffi.cdef('''int dist(const float* xyz, const int* pairs, float* distance_out,
                         float* displacement_out, const int n_frames, const int n_atoms,
                          const int n_pairs);''')
    ffi.cdef('''int angle(const float* xyz, const int* triplets, float* out,
                          const int n_frames, const int n_atoms, const int n_pairs);''')
    ffi.cdef('''int dihedral(const float* xyz, const int* quartets, float* out,
                             const int n_frames, const int n_atoms, const int n_pairs);''')

    here = os.path.dirname(os.path.abspath(__file__))
    libpath = find_library(here, 'geometry')
    print libpath
    if libpath is not None:
        C = ffi.dlopen(libpath)
        _HAVE_OPT = True
    else:
        _HAVE_OPT = False

if not _HAVE_OPT:
    warnings.warn('Optimized distance library was not imported sucessfully.')

import rg, internal, alignment
from .angle import *
from .distance import *
from .dihedral import *