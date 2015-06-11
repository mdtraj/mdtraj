MDTraj Programming Style
========================

PEP8
----

`PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ is the python style guide. Follow it. Having consistent style makes the code easy to read and easy to use. This is especially important for public methods (i.e. they should use ``underlines_to_separate_words`` and not ``camelCase``), but its also important for internal variable names, so that the code is easy to read and not jarring.

There are a few static code checkers that can help you with PEP8 compliance.
I like ``flake8``. There's also pylint. ::

  pip install flake8 pylint

There are a few stupid flake8 errors that you can overlook, but most of them should be corrected. The ones that are probably alright to skip are E501 (line too ling), E231, and a few others. Sometimes you just need more than 80 characters. Readability is what counts. Use your judgement. But err on the side of following PEP8. To ignore specific warnings with ``flake8``, you can run it as so. ::

  flake8 --ignore=errors=E501,E231 MDTraj/

Properties
----------

This is not java. Don't write ``get_myattribute`` and ``set_myattribute`` methods. Use the ``@property`` decorator. See the `astropy coding guidelines <http://docs.astropy.org/en/latest/development/codeguide.html#properties-vs-get-set>`_ for details.


Writing Tests
-------------
Untested code is broken code. So testing isn't really a feature of programming style, it's about correctness. Even tested code `might` be broken code, but untested code is by definition broken, and that's worse. Writing tests is annoying, but it pays off big time. `Nose <https://nose.readthedocs.org/en/latest/>`_ makes writing tests pretty easy, and MDTraj has some testing infrastructure to make it easier.

Tests should go in a file named ``test_modulename.py``, either in ``MDTraj/test/`` or in the ``test`` directory in a subpackage, like ``MDTraj/geometry/test``. If a subpackage only contains one test module, then it's fine to put the test module in the package directory directly, without making a ``test`` subdirectory. This is how the ``mdtraj.pdb`` subpackage's tests are organized currently.

MDTraj contains some infrastructure to help with testing. The most useful is the function ``mdtraj.testing.eq``,

.. autofunction:: mdtraj.testing.eq

.. autofunction:: mdtraj.testing.get_fn

Docstrings
----------

Every function should have a docstring. This is critical. Most of the information in the documentation is pulled directly from the docstrings, which means that the docstrings are really important. In order
to have sphinx properly parse the docstrings, they need to be written in the numpy format. A detailed description of the numpy docstring format is available `here <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Cython Docstrings
+++++++++++++++++
in cython code, the function signature cannot be introspected directly, due to limitations of python. We get around that by writing the signature manually on the first line of the docstring. For instance, the docstring for ``XTCTrajectoryFormat.read()`` is below. Note that the signature comes before the 1 line summary. ::

    """read(n_frames=None, stride=1, atom_indices=None)
    
    Read data from an XTC file

    Parameters
    ----------
    n_frames : int, None
        The number of frames you would like to read from the file.
        If None, all of the remaining frames will be loaded.
    stride : int, optional
        Read only every stride-th frame.
    atom_indices : array_like, optional
        If not none, then read only a subset of the atoms coordinates from the
        file. This may be slightly slower than the standard read because it required
        an extra copy, but will save memory.

    Returns
    -------
    xyz : np.ndarray, shape=(n_frames, n_atoms, 3), dtype=np.float32
        The cartesian coordinates, in nanometers
    time : np.ndarray, shape=(n_frames), dtype=np.float32
        The simulation time, in picoseconds, corresponding to each frame
    step : np.ndarray, shape=(n_frames), dtype=np.int32
        The step in the simulation corresponding to each frame
    box : np.ndarray, shape=(n_frames, 3, 3), dtype=np.float32
        The box vectors in each frame.
    """
 
Class Docstrings
++++++++++++++++
Make sure to include the constructor parameters in the class docstring, not the ``__init__`` docstring. This is more consistent with how the constructor is actually invoked, since users don't call ``cls.__init__(...)``, they call ``cls(...)``.

Also, it's a good idea to put an attributes section in the class docstring, even if the attributes are ``@property`` decorated functions.
