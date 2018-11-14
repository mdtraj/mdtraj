.. _installation::

Installation
============

We recommend that you install ``mdtraj`` with ``conda``. ::

  $ conda install -c conda-forge mdtraj

You can install also ``mdtraj`` with ``pip``, if you prefer. ::

  $ pip install mdtraj


Conda is a cross-platform package manager built especially for scientific
python. It will install ``mdtraj`` along with all dependencies from a
pre-compiled binary. If you don't have Python or the ``conda`` package
manager, we recommend starting with the `Anaconda Scientific Python
distribution <https://store.continuum.io/cshop/anaconda/>`_, which comes
pre-packaged with many of the core scientific python packages that MDTraj
uses (see below), or with the `Miniconda Python distribution
<http://conda.pydata.org/miniconda.html>`_, which is a bare-bones Python
installation.

MDTraj supports Python 2.7 or Python 3.4+ (recommended) on Mac, Linux, and
Windows.


Testing Your Installation
-------------------------

Running the tests is a great way to verify that everything is working. The test
suite uses `pytest <https://pytest.readthedocs.org/en/latest/>`_, which you can pick
up via ``pip`` if you don't already have it. ::

  pip install pytest

Then, to run the tests, execute the command ::

  pytesttests mdtraj -v

Compiling From Source
---------------------

To compile MDTraj from source, you'll need ``cython`` in addition to the
normal, runtime dependencies. Check ``devtools/conda-recipe/meta.yaml`` for
a complete list of build and run dependencies.

.. vim: tw=75
