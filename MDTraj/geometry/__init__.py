# This file is part of MDTraj.
#
# Copyright 2013 Stanford University
#
# MDTraj is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

__all__ = ['compute_distances', 'compute_angles', 'compute_dihedrals']

import os
import warnings
try:
    import cffi as _cffi
    from mdtraj.utils.ffi import find_library as _find_library
    _HAVE_OPT = None   # not sure if we have the library yet
except ImportError:
    warnings.warn('Optimized distance library requires the "cffi" package, '
                  'which is installable with easy_install or pip via '
                  '"pip install cffi" or "easy_install cffi".')
    _HAVE_OPT = False  # we definitely don't have the library

if _HAVE_OPT is not False:
    # lets try to open the library
    ffi = _cffi.FFI()
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
    ffi.cdef('''int kabsch_sander(const float* xyz, const int* nco_indices, const int* ca_indices,
                                  const int n_frames, const int n_atoms, const int n_residues,
                                  int* hbonds, float* henergies);''')

    _here = os.path.dirname(os.path.abspath(__file__))
    _libpath = _find_library(_here, 'geometry')
    if _libpath is not None:
        C = ffi.dlopen(_libpath)
        _HAVE_OPT = True
    else:
        _HAVE_OPT = False

if not _HAVE_OPT:
    warnings.warn('Optimized distance library was not imported sucessfully.')

import rg, internal, alignment, hbond
from .angle import *
from .distance import *
from .dihedral import *
