MDTraj
======

*Read, write and analyze MD trajectories with only a few lines of Python code.*

MDTraj is a python library that allows users to manipulate `molecular dynamics (MD) <http://en.wikipedia.org/wiki/Molecular_dynamics>`_ trajectories and perform a variety of analyses, including fast `RMSD <http://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions>`_,
`solvent accessible surface area <http://en.wikipedia.org/wiki/Accessible_surface_area>`_, hydrogen bonding, etc. A highlight of MDTraj is the wide variety of molecular dynamics trajectory file formats which are supported, including RCSB pdb, GROMACS xtc and trr, CHARMM / NAMD dcd, AMBER binpos, AMBER NetCDF, AMBER mdcrd, TINKER arc and :ref:`MDTraj HDF5 <HDF5FormatSpec>`.

MDTraj is based on numpy and is both easy to use and fast. Core routines like RMSD are written in C with explicit SSE vectorization and multicore parallelization. The RMSD code, in particular, is based on the Theobald QCP method and is 4x the speed of the original Theobald code and over 3x as fast as Theobald code modified to use GotoBLAS.

The library also ships with a flexible command-line application for converting
trajectories between formats. When you install MDTraj, the script will be
installed under the name ``mdconvert``.

.. raw:: html
 
  <div>
  <h2 style="display: inline; float:left; margin-left:5em"><a href="https://github.com/rmcgibbo/mdtraj/releases/latest">
  Download the Code</a></h2>
  <h2 style="display: inline; float:right; margin-right:5em"> <a href="examples/index.html"> See it in Action
  </a></h2>
  <div style="clear:both"></div>
  <div style="display:block; text-align:center;"><h2 style="display:inline;"><a href="https://github.com/rmcgibbo/mdtraj/issues">Get Involved</a></h2></div>
  </div>
  <br/>


Feedback
--------
The best way to report a bug or request a new feature is to make an issue
on github. Don't hesitate to fork the repository, make some changes, and 
submit a pull request!


--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
   :maxdepth: 1
   
   getting_started
   examples/index
   mdconvert
   whatsnew

API Reference
-------------

.. toctree:: 
   :maxdepth: 1
   
   api/load_functions
   api/trajectory_files
   api/classes
   api/analysis
   api/reporters
   api/utils


Developing
----------

.. toctree:: 
   :maxdepth: 1

   hdf5_format
   building_docs
   style


--------------------------------------------------------------------------------

License
-------
MDTraj is licensed under the Lesser GNU General Public License (LGPL v2.1+).
