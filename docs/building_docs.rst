.. _building-docs:

Building the documentation
--------------------------

The MDTraj documentation is built using `sphinx <http://sphinx-doc.org/>`_ and requires a few dependencies like `ipython <http://ipython.org/>`_ and `matplotlib <http://matplotlib.org/>`_ that you probably already have installed. We use `travis-ci <https://travis-ci.org/>`_ for continuous integration (running the tests), and also for building the documentation, which is built and pushed directly to Amazon S3 after every successful build.

Although `readthedocs <https://readthedocs.org/>`_ is a great tool, it doesn't have the flexibility we need for this project. We use sphinx's autodoc feature to generate documentation from docstrings, many of which are in compiled cython code. This means that building the documentation requires having a fully compiled version of MDTraj, which is not supported on readthedoc's servers. Furthermore, our documentation includes examples (with plots!) which are built directly with the documentation, and requires a fully functional MDTraj environment.

If you'd like to build the docs on your machine, you'll first need to install sphinx and numpydoc ::

    pip install sphinx numpydoc

You'll also need IPython, pandoc and runipy for the example notebooks ::

    sudo apt-get install pandoc
    pip install ipython runipy
  
Now, go back to the docs subdirectory in the main repository. The documentation will be built in the ``docs/_build`` subdirectory ::

    cd docs
    make html
