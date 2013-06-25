.. MDTraj documentation master file, created by
   sphinx-quickstart on Tue Jun 11 21:23:28 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MDTraj
======

MDTraj is a python library for loading, saving, and manipulating molecular
dynamics trajectory files.

Features
--------
MDTraj currently supports loading and saving trajectories in the following
formats.

- PDB
- GROMACS XTC and TRR
- CHARMM / NAMD DCD
- AMBER binpos, NetCDF
- MDTraj HDF5

The MDTraj geometry library gives you access to fast cartesian root-mean
square deviation (RMSD) distances, the calculation of bond lengths, bond
angles, dihedral angles, non-redundant internal coordinates, and radius
of gyration.

Source Code
-----------
The source code is licensed under the GPL, and available at: https://github.com/rmcgibbo/mdtraj

Feedback
--------
The best way to report a bug or request a new feature is to make an issue
on github. Don't hesitate to fork the repository, make some changes, and submit
a pull request!

API Reference
-------------

.. toctree::
   :maxdepth: 2

   getting_started
   file_formats