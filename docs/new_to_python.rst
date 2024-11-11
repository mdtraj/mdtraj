New to Python
=============

These instructions are intended for users from a biology, chemistry, or
physics perspective who are getting started with MDTraj and may be
unfamiliar with Python or programming.

Why Python?
-----------

MDtraj is not a standalone program. It's a software package that provides
functions for analyzing and working with molecular dynamics trajectories
in python; it's a set of building blocks.

Python is ideal because it gives you the ability to mix and match or
combine tools in new ways. You can truly develop analyses and
visualizations.  You can do all of this while building on top of a wealth
of high-quality open source packages for interactive scientific computing,
like `IPython <http://ipython.org/>`_, `numpy <http://www.numpy.org/>`_,
`scipy <http://scipy.org/>`_, `scikit-learn
<http://scikit-learn.org/stable/>`_, and `matplotlib
<http://matplotlib.org/>`_.


Installing Python
-----------------

There are multiple ways to install Python. We strongly recommend that you
use the `Anaconda distribution <http://continuum.io/downloads>`_, which is
aimed towards scientists. It comes pre-installed with most of the building
blocks you need for numerical and scientific computing, and with a great
package manager, ``conda``, for installing additional packages.

.. tip:: `Download <http://continuum.io/downloads>`_ Anaconda scientific
         python. We test against python 3.10, 3.11, and 3.12


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

  Python 3.5.1 :: Continuum Analytics, Inc.

If everything looks good, go ahead and install ``mdtraj`` using the ``conda``
package manager. ::

  conda install -c conda-forge mdtraj

The ``conda-forge`` "channel" provides a variety of packages supported by the
Python community.


Starting MDTraj with IPython/Jupyter
------------------------------------

If you're using Windows, open the Jupyter Notebook from the Start Menu. For
OS X or Linux, run the command ``jupyter notebook`` from a command line
prompt.  This should load up a new notebook in your browser.

From here, you can start by trying ``import mdtraj``, or copying a pasting
commands from our :ref:`examples <examples>` page, and getting familiar
with the environment. You're off to the races!


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

    `MDTraj Github Issue Tracker <https://github.com/mdtraj/mdtraj/issues>`_
        The development activity for MDTraj takes place openly on github. This
        is the best way to contact the developers directly, and participate
        in improving MDTraj for everyone.

.. vim: tw=75
