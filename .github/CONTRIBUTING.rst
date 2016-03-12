MDTraj Style Guide
==================

PEP8
----

Please follow `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_, the
python style guide.  There are a few static code checkers that can help you
with PEP8 compliance including "flake8". PyCharm includes pep8 checks and
can automatically fix pep8 problems.

When submitting a pull-request, don't reformat code unrelated to your work,
even if the other code is not pep8 compliant. This clutters the "diff" and
makes it difficult to review your changes.


Tests
-----

Writing unit tests in a dynamic language like Python is very important.
Writing tests is annoying, but it pays off big time. `Nose
<https://nose.readthedocs.org/en/latest/>`_ makes writing tests pretty
easy, and MDTraj has some testing infrastructure to make it easier. We will
probably refuse your pull-request if it doesn't include tests.

Tests should go in a file named ``test_modulename.py`` in
``mdtraj/tests/``.  Good tests each test a very specific piece of
functionality. To test multiple functionalities, write more tests. Unlike
normal programming, you can "copy-paste" code when writing tests.

Docstrings
----------

Every function should have a docstring. This is critical. Most of the
information in the documentation is pulled directly from the docstrings,
which means that the docstrings are really important. For the same reason,
it's really important that the formatting of docstrings is exactly correct.
In order for them to be parsed, they need to be written in the numpydoc
format. A detailed description of the numpy docstring format is available
`here
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Cython Docstrings
-----------------

in cython code, the function signature cannot be introspected directly, due
to limitations of python. We get around that by writing the signature
manually on the first line of the docstring. For instance, the docstring
for ``XTCTrajectoryFormat.read()`` is below. Note that the signature comes
before the 1 line summary. ::

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
 
.. vim: tw=75
