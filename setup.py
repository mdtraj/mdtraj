from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy

xtc = Extension('mdtraj.xtc',
    sources = ['MDTraj/xtc/src/xdrfile.c', 'MDTraj/xtc/src/trr2xtc.c',
               'MDTraj/xtc/src/xdrfile_trr.c', 'MDTraj/xtc/src/xdrfile_xtc.c',
               'MDTraj/xtc/xtc.pyx'],
    include_dirs = ["MDTraj/xtc/include/", 'MDTraj/xtc/', numpy.get_include()])
    
dcd = Extension('mdtraj.dcd',
    sources = ["MDTraj/dcd/src/dcdplugin.c", "MDTraj/dcd/dcd.pyx"],
    libraries=['m'],
    include_dirs = ["MDTraj/dcd/include/", 'MDTraj/dcd/', numpy.get_include()])


setup(name='MDTraj',
      packages=['mdtraj', 'mdtraj.pdb', 'mdtraj.xtc', 'mdtraj.dcd'],
      package_dir={'mdtraj':'MDTraj'},
      ext_modules=[xtc, dcd],
      cmdclass = {'build_ext': build_ext},
      include_package_data=True,
      package_data = {'mdtraj.pdb': ['data/*']})
