# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy

xtc = Extension('mdtraj.xtc',
    sources = ['MDTraj/xtc/src/xdrfile.c', 'MDTraj/xtc/src/xdrfile_xtc.c',
               'MDTraj/xtc/xtc.pyx'],
    include_dirs = ["MDTraj/xtc/include/", 'MDTraj/xtc/', numpy.get_include()])

trr = Extension('mdtraj.trr',
    sources = ['MDTraj/xtc/src/xdrfile.c', 'MDTraj/xtc/src/xdrfile_trr.c', 
               'MDTraj/xtc/trr.pyx'],
    include_dirs = ["MDTraj/xtc/include/", 'MDTraj/xtc/', numpy.get_include()])

dcd = Extension('mdtraj.dcd',
    sources = ["MDTraj/dcd/src/dcdplugin.c", "MDTraj/dcd/dcd.pyx"],
    libraries=['m'],
    include_dirs = ["MDTraj/dcd/include/", 'MDTraj/dcd/', numpy.get_include()])

binpos = Extension('mdtraj.binpos',
    sources = ['MDTraj/binpos/src/binposplugin.c', 'MDTraj/binpos/binpos.pyx'],
    include_dirs = ["MDTraj/binpos/include/", 'MDTraj/binpos/', numpy.get_include()])

setup(name='mdtraj',
      packages=['mdtraj', 'mdtraj.pdb', 'mdtraj.testing'],
      package_dir={'mdtraj':'MDTraj'},
      install_requires=['numpy', 'cython'],
      zip_safe=False,
      ext_modules=[xtc, trr, dcd, binpos],
      cmdclass = {'build_ext': build_ext},
      package_data = {'mdtraj.pdb': ['data/*'],
                      'mdtraj.testing': ["reference/*"]})
