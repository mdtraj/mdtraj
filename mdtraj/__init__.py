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

from mdtraj.formats.registry import _FormatRegistry
from mdtraj.formats.xtc import load_xtc
from mdtraj.formats.trr import load_trr
from mdtraj.formats.hdf5 import load_hdf5
from mdtraj.formats.lh5 import load_lh5
from mdtraj.formats.netcdf import load_netcdf
from mdtraj.formats.mdcrd import load_mdcrd
from mdtraj.formats.dcd import load_dcd
from mdtraj.formats.binpos import load_binpos
from mdtraj.formats.pdb import load_pdb
from mdtraj.formats.arc import load_arc
from mdtraj.formats.openmmxml import load_xml
from mdtraj.formats.prmtop import load_prmtop
from mdtraj.formats.psf import load_psf
from mdtraj.formats.mol2 import load_mol2
from mdtraj.formats.amberrst import load_restrt, load_ncrestrt
from mdtraj.formats.lammpstrj import load_lammpstrj
from mdtraj.formats.dtr import load_dtr
from mdtraj.formats.xyzfile import load_xyz

from mdtraj.core import element
from mdtraj._rmsd import rmsd
from mdtraj._lprmsd import lprmsd
from mdtraj.core.topology import Topology
from mdtraj.geometry import *
from mdtraj.core.trajectory import *
from mdtraj.nmr import *
import mdtraj.reporters

def test(label='full', verbose=2):
    """Run tests for mdtraj using nose.

    Parameters
    ----------
    label : {'fast', 'full'}
        Identifies the tests to run. The fast tests take about 10 seconds,
        and the full test suite takes about two minutes (as of this writing).
    verbose : int, optional
        Verbosity value for test outputs, in the range 1-10. Default is 2.
    """
    import mdtraj
    from mdtraj.testing.nosetester import MDTrajTester
    tester = MDTrajTester(mdtraj)
    return tester.test(label=label, verbose=verbose, extra_argv=('--exe',))
# prevent nose from discovering this function, or otherwise when its run
# the test suite in an infinite loop
test.__test__ = False


def capi():
    import os
    import sys
    module_path = sys.modules['mdtraj'].__path__[0]
    return {
        'lib_dir':  os.path.join(module_path, 'core', 'lib'),
        'include_dir': os.path.join(module_path, 'core', 'lib'),
    }
