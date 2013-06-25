Getting started
###############

Installing MDTraj
=================

MDTraj currently runs best on Python 2.7.x; earlier versions of Python are
not supported. We are working on support for Python 3.3+.

You can install MDTraj by cloning the repository from github::

	$ git clone git://github.com/rmcgibbo/mdtraj.git
	$ cd mdtraj
	$ python setup.py install
	
Or by installing it directly with `pip <http://www.pip-installer.org/>`_::

	pip install git+git://github.com/rmcgibbo/mdtraj.git
	
Dependencies
============

MDTraj requires a few external libraries to work properly. Some of them are
optional. You'll definitely need ``numpy`` and ``scipy``. You'll also need
a working C compiler and build environment, to compile the C-level extensions
that connect the python code to the bit-twiddling C libraries that read and write
most of these binary formats.

For the MDTraj HDF5 format, you'll need the python `HDF5 <http://www.pytables.org/moin>`_ library, ``pytables``.
For the AMBER NetCDF interface, you'll need the python `netcdf <https://pypi.python.org/pypi/netCDF4/1.0.4>`_ library,
``netcdf4-python``. Both of these are installed by default if you use the
`Enthought <https://www.enthought.com/products/epd/>`_ python distribution.
