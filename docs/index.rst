.. MDTraj documentation master file, created by
   sphinx-quickstart on Tue Jun 11 21:23:28 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MDTraj
======

MDTraj is a library for loading, saving, and manipulating molecular
dynamics trajectories in the python programming language.

Features
--------
MDTraj currently supports loading and saving trajectories in the RCSB PDB, GROMACS XTC and TRR, CHARMM / NAMD DCD, AMBER binpos, AMBER NetCDF, AMBER mdcrd, TINKER arc and :ref:`MDTraj HDF5 <HDF5FormatSpec>` formats.

MTraj is structured around a convent numpy-based ``Trajectory`` object.  Its RMSD (optimal cartesian root-mean-square deviation) library gives you access to `fast` parallel SSE3-based structural deviations using the quaternion characteristic polynomial (Theobald QCP) method at 4x the speed of the original Theobald code and over 3x as fast as Theobald code modified to use GotoBLAS.

The library also ships with a flexible command-line application for converting
trajectories between formats. When you install MDTraj, the script will be
installed under the name ``mdconvert``.

Documentation
-------------

.. toctree::
   :maxdepth: 1
   
   getting_started
   examples/index
   mdconvert
   hdf5_format
   file_formats
   loading
   trajectory
   topology
   rmsd
   io
   geometry
   compatibility
   reporters
   utils

Minimal Example
---------------
Start by loading up a trajectory from disk. MDTraj will automatically parse the file extension and use the appropriate loader. ::

  >>> import mdtraj as md
  >>> t = md.load('trajectory.xtc', top='trajectory.pdb')
  >>> print t
  <mdtraj.Trajectory with 100 frames, 22 atoms at 0x109f0a3d0>

To load files that don't contain topology information, like Gromacs XTC files,
we need to supply something with the ``top`` keyword argument that describes 
the topology, for example a PDB file.

If I'm interest in only a subset of the frames of the trajectory, I can slice it ::
  
  >>> # lets take a look at the first ten frames
  >>> print t[1:10]
  <mdtraj.Trajectory with 9 frames, 22 atoms at 0x109b91310>
  
  >>> # or maybe the last frame?
  >>> print t[-1]
  <mdtraj.Trajectory with 1 frames, 22 atoms at 0x109b97810>

There's a lot of information in the trajectory object. The most obvious is the
cartesian coordinates. They're stored as a numpy array under ``xyz``. All of
the distances in the :ref:`Trajectory` are stored in nanometers. The time unit
is picoseconds. Angles are stored in degrees (not radians). ::

  >>> print t.xyz.shape
  (100, 22, 3)
  >>> print np.mean(t.xyz)
  0.89365752249053032

  >>> # the simulation time (in picoseconds) of th first 10 frames
  >>> print t.time[0:10]
  array([ 0.002,  0.004,  0.006,  0.008,  0.01 ,  0.012,  0.014,  0.016,
        0.018,  0.02 ], dtype=float32)
  
  >>> # or the unitcell lengths in the last frame? (in nanometers of course)
  >>> t.unitcell_lengths[-1]
  array([ 2.,  2.,  2.], dtype=float32)
  

Saving the trajectory back to disk is easy. ::

  >>> # the hdf5 format stores the topology inside the file for convenience
  >>> t[::2].save('halftraj.h5')
  
  >>> # the format will be parsed based on the extension, or you can call the
  >>> # format-specific save methods
  >>> t[0:10].save_dcd('first-ten-frames.dcd')


The trajectory contains a reference to a topology object, which can come in handy. For example, if you want to save a copy of your trajectory with
only alpha carbons present, you can do that pretty easily. ::

  >>> atom_to_keep = [a.index for a in t.topology.atoms if a.name == 'CA']
  >>> t.restrict_atoms(atoms_to_keep)  # this acts inplace on the trajectory
  >>> t.save('CA-only.h5')

Source Code
-----------
The source code is licensed under the Lesser GNU General Pluglic License (v2.1+), and available at: https://github.com/rmcgibbo/mdtraj

Feedback
--------
The best way to report a bug or request a new feature is to make an issue
on github. Don't hesitate to fork the repository, make some changes, and 
submit a pull request!

Developer
---------

.. toctree::
   :maxdepth: 2

   building_docs
   style
