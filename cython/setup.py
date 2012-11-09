from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy

xtc = Extension('xtc',
    sources = ['xdrlib/xdrfile.c',
               'xdrlib/trr2xtc.c',
               'xdrlib/xdrfile_trr.c',
               'xdrlib/xdrfile_xtc.c',
               'xtc.pyx'],
    include_dirs = ["xdrlib/include/", numpy.get_include()])
    
setup(name='MDTraj',
      cmdclass = {'build_ext': build_ext},
      ext_modules = [xtc],
      include_package_data=True,
      package_data = {'mdtraj.pdb': ['data/*']})