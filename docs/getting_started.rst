.. _getting-started:

************
Installation
************

Platforms
=========

MDTraj runs on CPython 2.6, 2.7, and 3.2+. The developers generally use
CPython2.7 or 3.3, on both Mac and Linux platforms. Automated tests on Linux
and Windows are performed on every incremental update to the code, and release
tests are performed on Mac, Linux and Windows platforms.

Supported Hardware
------------------
Some components of MDTraj are written in C / C++ and make use of `SSE
intrinsics <http://en.wikipedia.org/wiki/Streaming_SIMD_Extensions>`_ for
maximum performance. The SSE4.1 instruction set (the latest we use) was released
in 2007, so most recent x86 machines should have no problem with this
requirement. Utilization of these features also requires using a suitably modern
compiler. This is only really a question on Linux, where gcc >= 4.3 is a
prerequisite.


Install with Conda
------------------
.. _install-with-conda:

`conda <http://www.continuum.io/blog/conda>`_ is a python package manager built for scientific python. Unlike ``easy_install`` or ``pip``, it handles binaries and binary dependencies, which are critical for most scientific workflows. If you're a ``conda`` user, you can install MDTraj by adding the omnia channel. If you're not a conda user, you should look into it. ::

To install MDTraj with conda, use the following commands ::

  $ conda config --add channels http://conda.binstar.org/omnia
  $ conda install mdtraj

.. note:: ``conda`` will automatically all of the tricky dependencies from binary packages automatically! This includes pytables / numpy / scipy! The easiest way to get conda is with the `Anaconda python distribution <https://store.continuum.io/cshop/anaconda/>`_.


Install with Pip
----------------

On a Mac or Linux machine, just run ::

  $ pip install mdtraj
  
Or, if you want the bleeding-edge source code, use ::

  $ pip install git+git://github.com/rmcgibbo/mdtraj.git

Install from Source
-------------------
Clone the source code repository from github::

  $ git clone git://github.com/rmcgibbo/mdtraj.git

If you don't have ``git``, you can download the source code as a zip file from
https://github.com/rmcgibbo/mdtraj/archive/master.zip. Then, in the directory containing the source code, you can install it with. ::

  $ python setup.py install

Dependencies
============

To use mdtraj, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Linux and Mac machines. Windows is not
        well supported.

    `Python <http://python.org>`_ >= 2.6
        The development package (``python-dev`` or ``python-devel``
        on most Linux distributions) is recommended.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        Numpy is the base package for numerical computing in python.

Optional packages:

    `SciPy <http://scipy.org>`_ >= 0.12.0
        We use scipy for loading and saving AMBER netcdf formatted
	trajectories. Note that due to a `serious bug
	<https://github.com/rmcgibbo/mdtraj/pull/466>` in the NetCDF
	support in scipy 0.11.0, a relatively modern version (0.12.0)
	or greater is required.

    `Pandas <http://pandas.pydata.org>`_ >= 0.9.0
        Some functionality, including mol2 parsing,  requires pandas.

    `cython <http://cython.org>`_ >= 0.19.0
        Compiling MDTraj from source requires cython, and a working
	C / C++ compiler.

    `PyTables <http://www.pytables.org/>`_ >= 2.4.0
        Working with HDF5 formatted trajectories requires the PyTables
	package.

Avoid Hassles with Anaconda or Canopy
-------------------------------------

The easiest way to get all of the dependencies is to install one of the 
pre-packaged scientific python distributes like `Enthought's Canopy 
<https://www.enthought.com/products/canopy/>`_ or `Continuum's Anaconda 
<https://store.continuum.io/>`_. These distributions already contain all of 
the dependences, and are distributed via 1-click installers for Windows, Mac 
and Linux.

.. note:: I (RTM) personally recommend Continuum's Anaconda. It's free, includes the **AWESOME** conda package manager, and quite simple to use.

On python26, two backported packages, ``argparse`` and ``importlib`` are required. You can install them with ``pip install argparse`` and ``pip install importlib``.

Manually Installing the Dependencies
------------------------------------

Linux
++++++
If you're on ubuntu and have root, you can install everything through your package manager. ::

    $ sudo apt-get install libhdf5-serial-dev python-dev python-numpy python-scipy python-nose python-setuptools cython python-numexpr python-tables python-pandas

Mac
+++
If you're on mac and want a package manager, you should be using `homebrew <http://mxcl.github.io/homebrew/>`_ and ``brews``'s python (see `this page <https://github.com/mxcl/homebrew/wiki/Homebrew-and-Python>`_ for details). The best way to install numpy and scipy with ``brew`` is:

  $ brew tap Homebrew/python
  $ brew install python
  $ brew install numpy
  $ brew install scipy
  $ brew install hdf5

Then, you can install the remaining packages with pip. ::

  $ pip install nose numexpr cython tables pandas

Windows
+++++++
Chris Gohlke maintains windows binary distributions for an ever-growing
set of python extensions on `his website <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.
Download and install the the installers for setuptools, nose, numpy, scipy, numexpr, pandas and tables.

Compiling Dependencies from source (no root needed)
---------------------------------------------------

If you don't already have a python installation you want to use, you can compile a new one. ::

  $ wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
  $ tar -xzvf Python-2.7.5.tgz
  $ cd Python-2.7.5
  $ ./configure --enable-shared --prefix=$HOME/local/python
  $ make
  $ make install

  $ export PATH=$HOME/local/python/bin:$PATH
  $ export LD_LIBRARY_PATH=$HOME/local/python/lib:$LD_LIBRARY_PATH

To  compile  the dependences  from  source,  you  need  to get  ``libhdf5``  and
``numpy``, which can  have some BLAS issues. I  recommend configuring everything
with  ``--prefix`` so  that you  don't get  your packages  mixed up  with system
packages. ::

  $ wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.11.tar.gz
  $ tar -xzvf hdf5-1.8.11.tar.gz
  $ cd hdf5-1.8.11
  $ ./configure --prefix=$HOME/opt/hdf5-1.8.11
  $ make
  $ make install

  $ export LD_LIBRARY_PATH=$HOME/opt/hdf5-1.8.11/lib:$LD_LIBRARY_PATH
  $ export PATH=$HOME/opt/hdf5-1.8.11/bin:$PATH

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
  $ pip install nose

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``pip`` if you don't already have it. ::

  pip install nose
  
Then, to run the tests, open a python shell and do ::

  >>> import mdtraj
  >>> mdtraj.test()

From the source directory, you can also run the tests with ``nosetests`` on
the command line
