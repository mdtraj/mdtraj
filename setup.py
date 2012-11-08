import glob
from setuptools import setup

setup(name='MDTraj',
      packages=['mdtraj', 'mdtraj.pdb', 'mdtraj.xtc', 'mdtraj.dcd'],
      include_package_data=True,
      package_data = {'mdtraj.pdb': ['data/*']})
