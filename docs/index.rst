MDTraj
======

*Read, write and analyze MD trajectories with only a few lines of Python code.*

MDTraj is a python library that allows users to manipulate `molecular dynamics (MD) <http://en.wikipedia.org/wiki/Molecular_dynamics>`_ trajectories.
Extensive :ref:`trajectory analysis <analysis>` routines are implemented. With MDTraj,
you can

 - Read and write from **every MD format imaginable** (pdb, xtc, trr, dcd, binpos,
   netcff, mdcrd, prmtop, ...)
 - Run **blazingly fast** RMSD calculations (4x the speed of the original
   `Theobald QCP <http://theobald.brandeis.edu/qcp/>`_).
 - Use tons of :ref:`analysis <analysis>` functions like bonds/angles/dihedrals,
   hydrogen bonding identification, secondary structure assignment, NMR observables.
 - **Lightweight API**, with a focus on **speed** and vectorized operations.

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
  <iframe id="player" type="text/html" width="500" height="300" style="display:block; margin:auto"
  src="http://www.youtube.com/embed/Lwy2Hdsr518"/>
  frameborder="0"></iframe>



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
   
   load_functions
   analysis
   api/trajectory_files
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
