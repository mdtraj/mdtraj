import glob
from setuptools import setup

setup(name='OpenMMPDB',
      packages=['OpenMMPDB'],
      include_package_data=True,
      package_data = {'OpenMMPDB': ['data/*']})
      #data_files=[('data', glob.glob('OpenMMPDB/data/*'))])
