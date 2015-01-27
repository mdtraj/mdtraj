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
import os
import mdtraj

from numpy.testing import Tester

class MDTrajTester(Tester):
    def _show_system_info(self):
        super(MDTrajTester, self)._show_system_info()
        print('mdtraj version %s' % mdtraj.version.version)
        print('mdtraj is installed in %s' % os.path.dirname(mdtraj.__file__))

        try:
            import tables
            print('tables version %s' % tables.__version__)
            print('tables hdf5 version %s' % tables.hdf5Version)
        except ImportError:
            print('tables is not installed')

        try:
            import netCDF4
            print('netCDF4 version %s' % netCDF4.__version__)
            print('netCDF4 lib version %s' % netCDF4.getlibversion())
        except ImportError:
            print('netCDF4 not installed')


