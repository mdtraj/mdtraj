from __future__ import print_function
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

        
