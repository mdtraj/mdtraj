***********
What's New?
***********

These are new features and improvements of note in each release.

v1.6 (Development)
------------------

MDTraj 1.6 contains a good mix of bug fixes and enhancements. Some
highlights include:

- Improved performance for ``compute_contacts`` (#995)
- Improved performance for ``Topology.select_pairs`` (#1000)
- Fast random access to xtc and trr files (#1038)
- xyz files support the ``__len__`` attribute (#998)
- ``segment_id`` is a new residue attribute (#1002)
- Expose ``FormatRegistry`` as a public api (#1039)
- Perform a heuristic check for valid unit cells when reading pdb files (#974)
- pdb file parsing uses the last model ``CONNECT`` records for bonds, not the first (#980)
- No longer force all warnings to be emitted (#1013 #1030)
- Always respect the ``force_overwrite`` argument in save methods (#878)
- Fix interop with ``scipy.cluster`` (#997)
- ``formats.hdf5.ensure_mode`` was removed (#990)


v1.5.1 (November 6, 2015)
-------------------------

MDTraj 1.5.1 is a small bugfix release to correct two issues introduced in the
immediately preceeding 1.5.0 release.

- A recent change (merged Nov 5) caused ``compute_chi4`` to compute chi3
  angles (#981).
- Revert changes in setup.py that resulted in a more confusing error when
  cython is not installed at build-time (#985).


v1.5 (November 6, 2015)
-----------------------

We're pleased to announce the 1.5 release of MDTraj. It contains new
features, improvements, and bug fixes. Highlights of the changes for this
version include:

- Faster histogramming method in RDF calculations when supported by numpy (#952)
- Improved support for mol2 reading (#945)
- Support for IPython/Jupyter 4 (#935)
- Improved support for Amber NetCDF writing (#939)
- Fix handling of periodic boundaries for distance calculations for general triclinic unit cells (#930)
- Support different reference and query indices for superposition and RMSD calculation (#915)
- Fix dcd reading bug under Windows (#905)
- Trajectories have a hash implementation (#898)
- Fixes for Hoomd (#900, #885)
- Support files (``devtools/``, ``setup.py``, ``.travis.yml``) are BSD licensed (#891, #893)
- Fixes for Lammpstrj (#861)
- Support for one letter amino acid codes (#871)
- Trajectory smoothing using a Buttersworth filter (#962)
- New functions for computing dihedral indices from a topology (#972)
- Improvements to build process (#955, #954, #941, #943, #942, #934)


v1.4.2 (June 9, 2015)
---------------------
- BUGFIX: Fix pytables inadvertently being moved to a required dependency


v1.4 (June 8, 2015)
-------------------
Version 1.4 is our best release yet! It contains many new features, performance improvements, and bug fixes.

Major highlights include:

- New function to calculate nematic order parameters (``compute_nematic_order``).
- Improved efficiency of generating RDF pairs.
- Add support for XYZ-format files.
- Fix parsing error with certain mol2 files.
- Support .pdb.gz files and make loading multiple pdb files more efficient.
- Fix use-after-free bug with DCD causing incorrect filenames.
- Update IPython-notebook trajectory viewer for IPython 3.0.
- Add support for the HOOMD-Blue XML topology format.
- Make virtual sites a new "element".
- Add 'NA' code to dssp for non-protein residues.
- Add support for CHARMM (Chamber) topologies in prmtop loader.
- Add methods to calculate more NMR J-couplings.
- Fix gro file unitcell handling.
- Enable .lammpstrj to parse custom column orders.
- Add read_as_traj method to all TrajectoryFile classes, making iterload work for all formats.

A total of 10 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.

Authors
~~~~~~~
* Kyle A. Beauchamp
* Anton Goloborodko +
* Matthew Harrigan
* Christoph Klein
* Robert T. McGibbon
* Tim Moore +
* Patrick Riley +
* Jason Swails
* Lee-Ping Wang
* Andrea Zonca +


v1.3 (February 25, 2015)
------------------------
- New functions to calculate various statistical mechanical properties
  (``unitcell_volumes``, ``dipole_moments``, ``static_dielectric``,
  ``isothermal_compressability_kappa_T``, ``thermal_expansion_alpha_P``,
  ``density``) (Kyle A. Beauchamp)
- Fix for PDB parser to handle more than 100K atoms. (Peter Eastman + ChayaSt)
- Include nitrogen atoms as h-bond acceptors in hydrogen bond detection (Gert Kiss)
- SSE4.1 support not required. The latest CPU feature now required is SSE3. (Robert T. McGibbon)
- New function to calculate radial distribution functions (``compute_rdf``) (Christoph Klein)
- Assorted bugfixes and improvements to documentation


v1.2 (December 1, 2014)
-----------------------
We're pleased to announce the 1.2 release of MDTraj! This release brings
minor changes to support the forthcoming release of MSMBuilder 3.

- Refactor RMSD code into a static library (Robert T. McGibbon)


v1.1 (November 10, 2014)
------------------------
We're pleased to announce the 1.1 release of MDTraj! This release brings
support for even more trajectory formats, and some new analysis features.

- New loader for CHARMM topology files: ``md.load_psf`` (Jason M. Swails)
- New loader for Desmond trajectory files (Teng Lin)
- New loader for Amber restart files (Jason M. Swails)
- New loader for Gromacs gro files (Robert T. McGibbon)
- New loader for LAMMPS trj files (Christoph Klein)
- New text-based :doc:`atom selection domain-specific language <atom_selection>`
  allowing natural querying of atoms as well as generation of equivalent
  python code for embedding in scripts or applications
  (Matthew P. Harrigan, Robert T. McGibbon)
- New ``md.compute_neighbors`` function to efficiently find nearby atoms (Robert T. McGibbon)
- ``md.shrake_rupley`` supports a new option to accumulate total SASA by residue
  (Robert T. McGibbon)
- Fix potential segmentation fault when reading corrupted XTC files.
  (Robert T. McGibbon)


v1.0.0 (September 7, 2014)
--------------------------
We're pleased to announce the 1.0 release of MDTraj! Our 1.0 release indicates
that MDTraj is stable enough to be used in production calculations, and that
we have a stronger commitment to backward compatibility. Two substantial new
features have been added since 0.9, but the API has remained quite stable.

- New interactive WebGl-based protein visualization in IPython notebook -- this
  feature is quite new and will continue to evolve throughout the 1.X release
  cycle.
- New ``md.compute_dssp`` function for DSSP secondary structure assignment.
- Multiple bugfixes in PDB parsing, including handling of ATOM serial's
  CONNECT records, support of .gziped files,
- Fix compilation errors on OSX and older linux platforms (gcc-4.1)
- ``compute_distances``, ``compute_angles``, ``compute_dihedrals`` now accept
  iterators for the indices argument.
- New ``Topology.select_atom_indices`` method.
- Ability to save b factors in PDB files.
- ``restrict_atoms`` has been deprecated, and replaced with ``atom_slice``.
- Better support for multi-chain proteins in dihedral methods.

Thanks to Robert T. McGibbon, Kyle A. Beauchamp, Lee-Ping Wang, Jason M. Swails,
ag1989, Carlos X. Hernandez, Matthew P. Harrigan and Christian Schwantes
for contributions.


v0.9.0 (June 10, 2014)
----------------------
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
