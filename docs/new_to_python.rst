Help! How do I get started?
===========================

These instructions are intended for users from a biology, chemistry, or
physics perspective who are getting started with MDTraj. It does not assume
prior experience with the Python programming language or the command line.


Why Python?
-----------

MDtraj is not a standalone program. It's a software package that provides
functions for analyzing and working with molecular dynamics trajectories
in python; **it's a set of building blocks.**

Python is ideal because it gives you, as a scientist, the ability to mix and
match, to combine tools in new ways, and to develop analyses and visualizations.
You can do all of this while building on top of a wealth of high-quality open
source packages for interactive scientific computing, like
`IPython <http://ipython.org/>`_, `numpy <http://www.numpy.org/>`_,
`scipy <http://scipy.org/>`_, `scikit-learn <http://scikit-learn.org/stable/>`_,
and `matplotlib <http://matplotlib.org/>`_.


Installing Python
-----------------
Because python is free, open source, and general purpose, there are a lot of
different ways to install it.

We strongly recommend that you use the `Anaconda scientific python
distribution <https://store.continuum.io/cshop/anaconda/>`_, because it comes
pre-installed with most of the building blocks you need for numerical and
scientific computing, and with a great package manager, ``conda``, for
installing additional packages. We also like the
`miniconda <http://conda.pydata.org/miniconda.html>`_ python distribution,
which is a little bit smaller and contains just python and ``conda``.

.. tip:: Here's the `download link <http://continuum.io/downloads>`_ for the Anaconda scientific python distribution.

When given the option to choose which version of python, we recommend the
latest, **python 3.4**.

On windows, we recommend using the 32-bit Anaconda installer, regardless of
whether or not your operating system is 32-bit or 64-bit. On OS X and Linux,
use the 64-bit python.


Installing MDTraj
-----------------
After installing python, you'll need to install ``mdtraj`` using the
command line.

- Windows: Find the ``Anaconda Command Prompt`` in the Start Menu
- OS X: Open the Terminal application in Applications / Utilities

From the command line, type ``python --version``. This should show the version
of python you have installed. For example, the command: ::

  $ python --version

prints out ::

  Python 3.4.2 :: Continuum Analytics, Inc.

If everything looks good, go ahead and install ``mdtraj`` using the ``conda``
package manager. ::

  conda install -c https://conda.binstar.org/omnia mdtraj


Starting MDTraj with IPython
----------------------------

If you're using Windows, open the IPython Notebook from the Start Menu. For
OS X or Linux, run the command ``ipython notebook`` from a command line prompt.
This should load up a new notebook in your browser.

From here, you can start by trying ``import mdtraj``, or copying a pasting
commands from our :ref:`examples <examples>` page, and getting familiar with
the environment. You're off to the races!


Resources and Getting Help
--------------------------

Here are some resources that you might find useful for getting started,
especially for scientists and scientific programming in Python.


Python resources:

    `Code Academy Introduction to Python <http://www.codecademy.com/en/tracks/python>`_
        This web course from Code Academy is a good introduction to basic
        programming in Python (without a particular focus on scientific
        works). It covers syntax, functions, control flow, and data structures.

    `Scientific Python Lecture Notes <https://scipy-lectures.github.io/>`_
        Tutorial material on the scientific Python ecosystem, a quick
        introduction to central tools and techniques. The different chapters
        each correspond to a 1 to 2 hours course with increasing level of
        expertise, from beginner to expert.

    `Lessons from the Software Carpentry <http://software-carpentry.org/lessons.html>`_
        The Software Carpentry is a non-profit organization whose members teach
        researchers basic software skills. All their teaching material is
        online and freely available. They might be hosting a workshop near you
        too!

Getting help:

    `MDTraj Discussion Forums <http://discourse.mdtraj.org/>`_
        You can ask questions in our discussion forums. It is generally quite
        responsive.

    `MDTraj Github Issue Tracker <https://github.com/mdtraj/mdtraj/issues>`_
        The development activity for MDTraj takes place openly on github. This
        is the best way to contact the developers directly, and participate
        in improving MDTraj for everyone.
