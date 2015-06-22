.. _installation::

************
Installation
************

We recommend that you install ``mdtraj`` with ``conda``. ::

  $ conda install -c omnia mdtraj

You can install also ``mdtraj`` with ``pip``, if you prefer. ::

  $ pip install mdtraj


Conda is a cross-platform package manager built especially for scientific python. It will install
``mdtraj`` along with all dependencies from a pre-compiled binary. If you don't have Python or the
``conda`` package manager, we recommend starting with the `Anaconda Scientific Python distribution
<https://store.continuum.io/cshop/anaconda/>`_, which comes pre-packaged with many of the core scientific
python packages that MDTraj uses (see below), or with the `Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`_, which is more bare-bones.

Supported Platforms
===================

Currently, we test and run ``mdtraj`` with Python 2.7, 3.3 and 3.4 on

- OS X 10.10 Yosemite
- x86-64 Ubuntu 12.04 LTS, 14.04 LTS
- x86-64 CentOS 5.x, 6.x
- 32-bit Python on 64-bit Windows Server 2012


Dependencies
============

To use ``mdtraj``, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Linux and OS X machines. We make our best
        effort to support Windows where possible.

    `Python <http://python.org>`_ >= 2.7
        The development package (``python-dev`` or ``python-devel``
        on most Linux distributions) is recommended.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        Numpy is the base package for numerical computing in python.

Optional packages:

    `SciPy <http://scipy.org>`_ >= 0.12.0
        We use scipy for loading and saving AMBER netcdf formatted
        trajectories.

    `Pandas <http://pandas.pydata.org>`_ >= 0.9.0
        Some functionality, including mol2 parsing,  requires pandas.

    `PyTables <http://www.pytables.org/>`_ >= 2.4.0
        Working with HDF5 formatted trajectories requires the PyTables
        package.


Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``pip`` if you don't already have it. ::

  pip install nose

Then, to run the tests, execute the command ::

  nosetests mdtraj -v

Compiling From Source
=====================

To compile MDTraj from source, you'll need ``cython``, ``numpy``, and an appropriate
compiler toolchain for your platform. For linux and OS X, we support the gcc and
clang compilers. For Windows, we support Microsoft Visual Studio 2008 and
Visual Studio 2010. The mingw compiler toolchain is unsupported.
