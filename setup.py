from setuptools import setup, Extension

xtc = Extension('mdtraj._xtc',
                sources = ['mdtraj/xtc/src/xdrfile.c',
                           'mdtraj/xtc/src/trr2xtc.c',
                           'mdtraj/xtc/src/xdrfile_trr.c',
                           'mdtraj/xtc/src/xdrfile_xtc.c'],
                include_dirs = ["mdtraj/xtc/src/include/"])

# dcd reader
dcd = Extension('mdtraj._dcd',
                sources = ["mdtraj/dcd/src/dcdplugin_s.c"],
                libraries=['m'],
                include_dirs = ["mdtraj/dcd/src/include/",
                                "mdtraj/dcd/src/"])
setup(name='MDTraj',
      packages=['mdtraj', 'mdtraj.pdb', 'mdtraj.xtc', 'mdtraj.dcd'],
      ext_modules=[xtc, dcd],
      include_package_data=True,
      package_data = {'mdtraj.pdb': ['data/*']})
