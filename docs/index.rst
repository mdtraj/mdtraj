MDTraj
======

Read, write and analyze MD trajectories with only a few lines of Python
code.

MDTraj is a python library that allows users to manipulate `molecular
dynamics (MD) <http://en.wikipedia.org/wiki/Molecular_dynamics>`_
trajectories.  Features include:

 - Wide MD format support, including ``pdb``, ``xtc``,
   ``trr``, ``dcd``, ``binpos``, ``netcdf``, ``mdcrd``, ``prmtop``, and
   more.
 - Extremely fast RMSD calculations (4x the speed of the original `Theobald
   QCP <http://theobald.brandeis.edu/qcp/>`_).
 - Extensive :ref:`analysis <analysis>` functions including those that
   compute bonds, angles, dihedrals, hydrogen bonds, secondary structure,
   and NMR observables.
 - Lightweight, Pythonic API.

MDTraj includes a command-line application, ``mdconvert``, for converting
trajectories between formats.


.. toctree::
   :maxdepth: 1

   installation
   new_to_python
   examples/index
   whatsnew

API Reference
-------------

.. toctree::
   :maxdepth: 1

   load_functions
   analysis
   atom_selection
   api/trajectory_files
   api/reporters
   api/utils
   mdconvert


Developing
----------

.. toctree::
   :maxdepth: 1

   hdf5_format



Citation |DOI for Citing MDTraj|
--------------------------------

MDTraj is research software. If you make use of MDTraj in scientific
publications, please cite it. The BibTeX reference is

::

    @article{McGibbon2015MDTraj,
        title = {MDTraj: A Modern Open Library for the Analysis of
        Molecular Dynamics Trajectories},
        author = {McGibbon, Robert T. and Beauchamp, Kyle A. and Harrigan,
        Matthew P. and Klein, Christoph and Swails, Jason M. and
        Hern{\'a}ndez, Carlos X.  and Schwantes, Christian R. and Wang,
        Lee-Ping and Lane, Thomas J. and Pande, Vijay S.},
        journal = {Biophysical Journal},
        volume = {109},
        number = {8},
        pages = {1528 -- 1532},
        year = {2015},
        doi = {10.1016/j.bpj.2015.08.015}
    }


License
-------

MDTraj is licensed under the Lesser GNU General Public License (LGPL
v2.1+).



.. |DOI for Citing MDTraj| image:: https://img.shields.io/badge/DOI-10.1016%2Fj.bpj.2015.08.015-blue.svg
   :target: http://doi.org/10.1016/j.bpj.2015.08.015

.. vim: tw=75
