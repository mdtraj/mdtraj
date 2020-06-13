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


"""The mdtraj package contains tools for loading and saving molecular dynamics
trajectories in a variety of formats, including Gromacs XTC & TRR, CHARMM/NAMD
DCD, AMBER BINPOS, PDB, and HDF5.
"""
import numpy as _  # silence cython related numpy warnings, see github.com/numpy/numpy/pull/432

from .formats.registry import FormatRegistry
from .formats.xtc import load_xtc
from .formats.trr import load_trr
from .formats.hdf5 import load_hdf5
from .formats.lh5 import load_lh5
from .formats.netcdf import load_netcdf
from .formats.mdcrd import load_mdcrd
from .formats.dcd import load_dcd
from .formats.binpos import load_binpos
from .formats.pdb import load_pdb
from .formats.arc import load_arc
from .formats.openmmxml import load_xml
from .formats.prmtop import load_prmtop
from .formats.psf import load_psf
from .formats.mol2 import load_mol2
from .formats.amberrst import load_restrt, load_ncrestrt
from .formats.lammpstrj import load_lammpstrj
from .formats.dtr import load_dtr, load_stk
from .formats.xyzfile import load_xyz
from .formats.hoomdxml import load_hoomdxml
from .formats.tng import load_tng


from .core import element
from ._rmsd import rmsd, rmsf
from ._lprmsd import lprmsd
from .core.topology import Topology, Single, Double, Triple, Amide, Aromatic
from .geometry import *
from .core.trajectory import *
from .nmr import *

from . import reporters

def capi():
    import os
    import sys
    module_path = sys.modules['mdtraj'].__path__[0]
    return {
        'lib_dir':  os.path.join(module_path, 'core', 'lib'),
        'include_dir': os.path.join(module_path, 'core', 'lib'),
    }
