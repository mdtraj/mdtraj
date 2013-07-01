.. _getting-started:

Getting started
###############

Installing MDTraj
=================

MDTraj currently runs best on Python 2.7.x; earlier versions of Python are not
supported. We are working on support for Python 3.3+. MDTraj is developed and
tested on mac and linux platforms. It probably works on windows, but might
require some hacking.

Easy Way
--------

Just run ::

  pip install mdtraj

Medium Way
----------
Clone the source code repository from github::

  $ git clone git://github.com/rmcgibbo/mdtraj.git

If you don't have ``git``, you can download the source code as a zip file from
https://github.com/rmcgibbo/mdtraj/archive/master.zip. Then, in the directory containing the source code, you can install it with. ::

  $ python setup.py install

Running the tests
=================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``pip`` if you don't already have it. ::

  pip install nose
  
Then, to run the tests, open a python shell and do ::

  >>> import mdtraj
  >>> mdtraj.test()

From the source directory, you can also run the tests with ``nosetests`` on
the command line

Dependencies
============

MDTraj requires a few external libraries to work properly. Some of them are
optional. You'll definitely need ``numpy`` and ``scipy``. You'll also need
a working C compiler and build environment, to compile the C-level extensions
that connect the python code to the bit-twiddling C libraries that read and 
write most of these binary formats.

Easy Way
--------

The easiest way to get all of the dependencies is to install one of the 
pre-packaged scientific python distributes like `Enthought's Canopy 
<https://www.enthought.com/products/canopy/>`_ or `Continuum's Anaconda 
<https://store.continuum.io/>`_. These distributions already contain all of 
the dependences, and are distributed via 1-click installers for Windows, Mac 
and Linux. You won't have to install anything except mdtraj, which can be done, as shown above, with ::

  pip install git+git://github.com/rmcgibbo/mdtraj.git

Medium Way
----------

Linux
++++++
If you're on ubuntu and have root, you can install everything through your package manager. ::

  $ sudo apt-get install libhdf5-serial-dev python-dev python-numpy python-nose python-setuptools cython python-numexpr python-tables netcdf-bin libnetcdf-dev python-netcdf python-networkx python-netcdf

Mac
+++
If you're on mac and want a package manager, you should be using `homebrew <http://mxcl.github.io/homebrew/>`_ and ``brews``'s python (see `this page <https://github.com/mxcl/homebrew/wiki/Homebrew-and-Python>`_ for details). The best way to install numpy and scipy with ``brew`` is to use
samueljohn's tap. ::

  $ brew tap samueljohn/python
  $ brew install python
  $ brew install numpy
  $ brew install scipy
  $ brew install hdf5
  $ brew install netcdf

Then, you can install the remaining packages with pip. ::

  $ pip install nose numexpr cython tables
  $ pip install netcdf4

Harder Way : Compiling from source (no root needed)
---------------------------------------------------

If you don't already have a python installation you want to use, you can compile a new one. ::

  $ wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
  $ tar -xzvf Python-2.7.5.tgz
  $ cd Python-2.7.5
  $ ./configure --prefix=$HOME/local/python
  $ make
  $ make install

  $ export PATH=$HOME/local/python/bin:$PATH

To compile the dependences from source, you need to get ``libhdf5``, ``libnetcdf``, and ``numpy``, which can have some BLAS issues. I recommend configuring everything with ``--prefix`` so that you don't get your packages mixed up with system packages. NetCDF4 depends on hdf5, so you need to install HDF5 first. ::

  $ wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.11.tar.gz
  $ tar -xzvf hdf5-1.8.11.tar.gz
  $ cd hdf5-1.8.11
  $ ./configure --prefix=$HOME/opt/hdf5-1.8.11
  $ make
  $ make install

  $ export LD_LIBRARY_PATH=$HOME/opt/hdf5-1.8.11/lib:$LD_LIBRARY_PATH
  $ export PATH=$HOME/opt/hdf5-1.8.11/bin:$PATH

  $ cd ..
  $ wget http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-4.3.0.tar.gz
  $ tar -xzvf netcdf-4.3.0.tar.gz
  $ cd netcdf-4.3.0
  $ CFLAGS="-I$HOME/opt/hdf5-1.8.11/include -L$HOME/opt/hdf5-1.8.11/lib"
  $ ./configure --prefix=$HOME/opt/netcdf-4.3.0
  $ make
  $ make install

  $ export LD_LIBRARY_PATH=$HOME/opt/netcdf-4.3.0/lib:$LD_LIBRARY_PATH
  $ export PATH=$HOME/opt/netcdf-4.3.0/bin:$PATH

You'll probably want to add those ``export`` statements to your bashrc too.

If you don't have ``easy_install`` or ``pip`` yet, you can get them with ::

  $ wget http://pypi.python.org/packages/source/s/setuptools/setuptools-0.6c11.tar.gz
  $ tar -xzvf setuptools-0.6c11.tar.gz
  $ cd setuptools-0.6c11.tar.gz
  $ python setup.py install
  $ easy_install pip

Now you're home free ::

  $ pip install numpy
  $ pip install scipy
  $ pip install cython
  $ pip install numexpr
  $ pip install tables
  $ pip install netcdf4
  $ pip install nose
