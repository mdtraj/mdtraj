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

from __future__ import print_function, division
from mdtraj.registry import _FormatRegistry
from mdtraj.xtc import XTCTrajectoryFile, load_xtc
from mdtraj.trr import TRRTrajectoryFile, load_trr
from mdtraj.hdf5 import HDF5TrajectoryFile, load_hdf5
from mdtraj.lh5 import LH5TrajectoryFile, load_lh5
from mdtraj.netcdf import NetCDFTrajectoryFile, load_netcdf
from mdtraj.mdcrd import MDCRDTrajectoryFile, load_mdcrd
from mdtraj.dcd import DCDTrajectoryFile, load_dcd
from mdtraj.binpos import BINPOSTrajectoryFile, load_binpos
from mdtraj.pdb import PDBTrajectoryFile, load_pdb
from mdtraj.arc import ArcTrajectoryFile, load_arc
from mdtraj.openmmxml import load_xml


from mdtraj._rmsd import rmsd
from mdtraj.topology import Topology
from mdtraj.geometry import *
from mdtraj.trajectory import *

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
