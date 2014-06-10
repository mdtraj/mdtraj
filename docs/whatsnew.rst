***********
What's New?
***********

These are new features and improvements of note in each release.

v0.9.0 (June 10, 2014)
-----------------------
- Brand new ``nmr`` library that includes transparent python interfaces to
  SHIFTX2, PPM and SPARTA+ for chemical shifts, as a library for scalar
  couplings (J) using the Karplus relation.
- New ``lprmsd`` distance metric for linear-programming root mean square
  deviations which optimizes over the label correspondence between
  exchangeable atoms in the two conformations.
- New ``wernet_nilsson`` function for hydrogen bond identification.
- New parser for ``mol2`` format files.
- Many new convenience methods on ``md.Topology``, including ``to_bondgraph``
  to create a NetworkX graph from a topology.
- New ``compute_drid`` function for calculation of distribution of
  reciprocal inter-atomic distances (DRID) distance metric
- Core geometry routines ``compute_angles`` and ``compute_dihedrals`` now
  respect periodic boundary conditions via a substantial internal refactoring
  of the geometry library. They also have significantly improved numerical
  stability.
- Numerous bugfixes, including fixing potential segfaults with ``md.rmsd`` and
  the NetCDF parser as well as increased compliance for AMBER .prmtop and
  TINKER .arc parsers.
- Many internal changes to hardware detection code, ensuring that compiled
  binaries run appropriately on any platform, including those that don't support
  modern CPU features like SSE4.
- Major improvements to our automated testing framework. Every pull request
  and commit to MDTraj is now being tested across a matrix of 4 different
  python versions on linux as well as python3 on Windows.
- A number of brand new example IPython notebooks on the website demonstrating
  all of these new features!


v0.8.0 (March 10, 2014)
-----------------------
- New parser for AMBER PRMTOP topology files.
- Removed dependency on netCDF4 and the c libnetcdf. We're now exclusively using
  the pure python NetCDF3 implementation in ``scipy.io``, which is now a dependency.
- Removed dependency on ``simtk.unit`` as an external package
- Fixed a behavior where "default" unit cell dimensions were being saved in
  trajectories without periodic boundary conditions in XTC, DCD and TRR, which
  when loaded up later were interpreted as being "real" periodic boundary conditions.
- Better ResSeq preservation in HDF5 files.
- More detailed ``repr`` and ``str`` on ``Trajectory``.
- Load pdb files directly from a URL.
- Unicode fixes for python3.
- Bugfixes in OpenMM reporters
- New theme for the documentation with IPython notebooks for the examples
- Improvements to ``DCD seek()``
- Reorganized the internal layout of the code for easier navigation, IPython
  tab completion.

Thanks to everyone who contributed to this release: Robert T. McGibbon,
Kyle A. Beauchamp, Carlos Hernandez, TJ Lane, Gert Kiss, and Matt Harrigan.

v0.7.0 (February 21, 2014)
--------------------------
- New geometry functions ``md.compute_contacts`` and ``md.geometry.squareform`` for residue-residue contact maps
- Fix segfault in ``md.rmsd`` when using the optional ``atom_indices`` kwarg
- ``md.compute_phi``, ``md.compute_psi``, and ``md.compute_omega`` now return the correct atom indices, as their docstring always said.
- Topology ``Element`` instances are now properly immutable
- Small bugfixes to ``baker_hubbard``, and better docstring
- Automatic installation of ``pandas`` and ``simtk.unit`` via setuptools' ``install_requires``.
- Small bugfix to mdcrd loading with stride
- ``superpose`` now correctly translates the final structure, and doesn't recenter the reference structure

v0.6.1 (February 11, 2014)
--------------------------
- ``Trajectory.join(discard_overlapping_frames=True)`` is criterion for detecting overlapping frames is more realistic
- We now support installation via conda, and are supplying conda binaries
- ``md.load()`` is much faster when loading multiple trajectory files
- Bug-fixes for pandas 0.13.0 release, detection of zinc atoms in PDB files
- Geometry functions are more resilient to segfaults from bad user parameters
- Fix intermittent RMSD segfaults from invalid memory access
- Fix RMSD centering bug with memory alignment after restrict_atoms

v0.6.0 (January 21, 2014)
-------------------------
- ``md.rmsd()`` signature changed to be more understandable
- All file objects now have a ``__len__`` function.
- Small bugfixes related to vsites.

v0.5.1 (January 4, 2014)
------------------------
- Minor bug fix when no dihedrals match specification
- Add ``__str__`` to Topology parts
- More examples sections in docstrings

v0.5.0 (January 3, 2014)
------------------------
- Numerous bug fixes
- Much improved coverage of the test suite.
- Removed cffi dependency for accelerated geometry code
- Faster multi-trajectory loading
- MSMBuilder2 LH5 format support
- Change license from GPL to LGPL
- More convenience methods on Topology
- PDB writer writes connect records
- Hydrogen bond identification with ``baker_hubbard``
- Rotation/translation to superpose trajectories
- New RMSD API. It's much simpler and much more memory efficient
- Full support for computing all of the chi angles
- Add seek/tell methods to all of the trajectory file objects
- New top level memory efficient ``iterload`` method for chunked trajectory loading
