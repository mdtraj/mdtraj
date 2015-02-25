MDTraj
======

*Read, write and analyze MD trajectories with only a few lines of Python code.*

MDTraj is a python library that allows users to manipulate `molecular dynamics (MD) <http://en.wikipedia.org/wiki/Molecular_dynamics>`_ trajectories.
Extensive :ref:`trajectory analysis <analysis>` routines are implemented. With MDTraj,
you can

 - Read and write from **every MD format imaginable** (``pdb``, ``xtc``, ``trr``,
   ``dcd``, ``binpos``, ``netcdf``, ``mdcrd``, ``prmtop``, ...)
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
      <h2 style="display: inline; float:left; margin-left:5em">
          <a href="https://github.com/mdtraj/mdtraj/releases/latest">
          Download the Code</a>
      </h2>
      <h2 style="display: inline; float:right; margin-right:5em">
          <a href="examples/index.html">
          See it in Action</a>
      </h2>
      <div style="clear:both"></div>
      <h2 style="display: inline; float:left; margin-left:7em">
          <a href="http://discourse.mdtraj.org/">
          Get Help</a>
      </h2>
      <h2 style="display: inline; float:right; margin-right:7em">
          <a href="https://github.com/mdtraj/mdtraj/issues">
          Get Involved</a>
      </h2>
  </div>
  <br/>
  <iframe id="player" type="text/html" width="500" height="300" style="display:block; margin:auto"
  src="http://www.youtube.com/embed/Lwy2Hdsr518"/>
  frameborder="0"></iframe>
  <br/><br/>


.. raw:: html

   <div style="display:none">


--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
   :maxdepth: 1

   installation
   new_to_python
   examples/index
   whatsnew
   faq
   Discussion Forums <http://discourse.mdtraj.org>

API Reference
-------------

.. toctree::
   :maxdepth: 1

   load_functions
   atom_selection
   analysis
   api/trajectory_files
   api/reporters
   api/utils
   mdconvert


Developing
----------

.. toctree::
   :maxdepth: 1

   hdf5_format
   building_docs
   style


--------------------------------------------------------------------------------

.. raw:: html

   </div>


Citation |DOI for Citing MDTraj|
--------------------------------

MDTraj is research software. If you make use of MDTraj in scientific publications, please cite it. The BibTeX reference is

::

    @article {
        author = {McGibbon, Robert T. and Beauchamp, Kyle A. and Schwantes, Christian R. and Wang, Lee-Ping and Hern{\'a}ndez, Carlos X. and Harrigan, Matthew P. and Lane, Thomas J. and Swails, Jason M. and Pande, Vijay S.},
        title = {MDTraj: a modern, open library for the analysis of molecular dynamics trajectories},
        year = {2014},
        doi = {10.1101/008896},
        publisher = {Cold Spring Harbor Labs Journals},
        journal = {bioRxiv}
    }


License
-------
MDTraj is licensed under the Lesser GNU General Public License (LGPL v2.1+).



.. |DOI for Citing MDTraj| image:: https://img.shields.io/badge/DOI-10.1101%2F008896-blue.svg
   :target: http://doi.org/10.1101/008896
