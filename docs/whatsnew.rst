What's New?
===========

These are new features and improvements of note in each release.

v1.10.2 (December 11, 2024)
--------------------------

- Add netCDF4 as a dependency in ``setup.py`` (#1944)
- Support NumPy 2 (#1928)

Authors
~~~~~~~

- Jeremy Leung
- Jessica A. Nash

v1.10.1 (October 17, 2024)
--------------------------

- Update RMSD module for Cython 3.0 (#1917)
- Fix edge case with float rounding in angle computation (#1841)
- Add atom selection via atom_indices to SASA computation (#1844)
- ReadTheDocs migration (#1888)
- Fix netCDF import (#1894)
- Use standard library `ast.unparse` instead of third-party `astunparse` (#1898)
- Allow appending to HDF5 files (#1899)
- Remove deprecated `oldest-supported-numpy` (#1902)
- Remove remaining uses of `six` (#1908)
- Replace `xdrlib` with `mda_xdrlib` (#1911)
- Add `mda-xdrlib` to `setup.py` (#1913)
- Set up `versioneer` (#1909)
- Remove `distutils` from tests (#1916)
- Remove unused `xdrlib` import from Python layer (#1919)
- Fix PDB "vmd convention" false-positive (#1904)

A total of 11 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.

- Charlie Laughton
- Guillermo Pérez-Hernández
- Sukrit Singh
- Jeremy Leung
- Sander Roet
- Chris Jones+
- David W.H. Swenson
- Matthew W. Thompson
- Robert Schütz+
- Patrick Kunzmann
- Augustin Zidek+

v1.10.0 (May 30, 2024)
----------------------

- Remove TNG support (#1875)
- Drop support for long-unsupported upstreams (#1838)
- Drop `distutils` (#1834)
- Lint codebase (#1874)
- Allow `mdtraj.load()`` to read various non-standard PDB files (#1849)
- Use netCDF4 as default in `NetCDFTrajectoryFile`, fallback to SciPy implementation (#1878)
- Add dependabot for updating actions (#1865)

A total of 4 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.

- Jinzhe Zeng+
- Matthew W. Thompson
- Jeremy Leung
- Jessica A. Nash+

v1.9.9 (July 22, 2023)
---------------------

- Pin `Cython~=0.29` in `pyproject.toml` (#1801)
- Remove "stale" bot (#1758)
- Announce maintainer update (#1795)
- Use `oldest-supported-numpy` in `pyproject.toml` (#1751)

A total of 4 people contributed to this release, including three new contributors.
People with a "+" by their names contributed a patch for the first time.

In-Ho Yi+
Sukrit Singh+
Jeremy Leung+
Matthew W. Thompson

v1.9.8 (July 2, 2023)
---------------------

- Patch util_arm.h for M1 (#1694)
- Avoid side effects in mdtraj.load (#1706)
- Add PDB chainID support (#1715)
- Remove imports from ``simtk`` namespace (#1698)
- Implement reading and writing PDBx/mmCIF (#1718)
- Fix `compute_inertia_tensor` docstring (#1721)
- Fix typo in docs (#1725)
- Adds lower-level `compute_distances_core` function (#1728)
- Implement generic C code paths for non-accelerated architectures (#1727)
- Add note about project maintinance status (#1738)
- Force new residue when creating new chain (#1740)
- Warn if ``atom_indices`` are not monotonically increasing (#1766)
- Delete notes in issue-template by (#1778)
- Fix hashing when unitcell data is not present (#1781)
- Test on Python 3.11 (#1758)
- Constrain GSD to version 2 (#1790)

A total of 15 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.

- Robert T. McGibbon
- Sander Roet
- Peter Eastman
- Charlie Laughton
- Christoph Klein
- Luca Naef +
- Jacob Lumpkins +
- Ray A. Matsumoto
- Matthew W. Thompson
- Moritz Hoffmann +
- German P. Barletta +
- Masataka Yamauchi +
- Samuel Lotz +
- Toni G +
- Bojun Liu +

v1.9.7 (November 12, 2021)
-----------------------

 - Replace discourse link with gitter link (904338b799842c103dcb9e306e6878a739a4d39f)
 - Faster load function and more homogeneus file parser interface (#1648)
 - Reduce memory usage of ``rdf_t`` (#1661)
 - Update to new OpenMM namespace (#1668)
 - Fix ``compute_contact`` bug with glycine sidechains (#1674)
 - Fix errors in RMSF documentation (#1676)
 - Handle path-like objects in place of filenames (#1680)
 - Pin NumPy in pyproject.toml for Python 3.10 (#1681)
 - Fix compilation errors on M1 Macs (#1684)
 - Update Python versions in CI (#1686)
 - Fix jupyter runners in CI (#1687)

A total of 10 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.

- Robert T. McGibbon
- Sander Roet
- Charlie Laughton +
- Alexander Matthew Payne +
- Tucker Burgin +
- Maurice Karrenbrock
- Luca Naef +
- Jacob Lumpkins +
- Ray A. Matsumoto
- Matthew W. Thompson


v1.9.6 (April 20, 2021)
-----------------------
 - Fix compatibility with all versions of astunparse (03753d736e30f15f8f210434e689e9ff664bb611)
 - Rework CI to be simpler and more maintainable
 - Fix deployment of the documentation to the website
 - Don't use serials for more than 1 chain (#1612)
 - Added ``enforcePeriodicBox`` option for HDF5Reporter (#1622)
 - Add time-dependent distance and RDF functions (#1633)
 - Add ``select`` option to ``compute_center_of_mass()`` (#1640)
 - ``Topology.join()`` can optionally updates resSeq of new residues (#1639)

A total of 7 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.

- Robert T. McGibbon
- Sander Roet
- Ray A. Matsumoto
- Maurice Karrenbrock +
- Kirill Shmilovich
- Heyi Liang +
- Matthew W. Thompson


v1.9.5 (Jan 4, 2021)
--------------------

 - Fix memory allocation when opening XTC trajectory on OSX Catalina (#1594)
 - Write out serial instead of index in pdbs (#1584)
 - Fix residue idx sliced traj (#1586)
 - Update shift_wrappers.py (#1579)
 - Rsmd atom_indices checks fix (#1571)
 - Port to aarch64 (#1562)
 - Add compatibility with pandas 1.0


v1.9.4 (May 10, 2020)
-----------------------

- Update some pandas calls for v1.0 (#1536)
- Fix TRR file offset (#1534)
- Update selection for Python 3.8 compatibility (#1523)
- Ensure bonds exist before using them (#1512, #1513)
- Let compute_displacements handle empty atom_pairs (#1515)
- Add GSD reader and writer (#1494)
- Fix stride parameter for .netcdf files (#1501)
- Ensure that the license file is packaged in the sdist (#1498)
- Right-justify atom symbol when writing PDB files (#1459)
- Add calculations for shape metrics (#1471)
- Fix residue parsing in MOL2 reader (#1490)
- Set up "stale" bot
- Use AZP for CI (#1484, #1536)
- Fix leaving malformed TRR files open (#1482)
- Fix various OpenMP issues (#1476, #1477, #1488, #1508, #1529)
- Add gyration tensor calculation (#1467)
- Fix some type conversions (#1466, #1511)
- Remove bundled dependencies astor and pyparsing (#1452)
- Correct ordering in hoomdxml files (#1453)

Authors
~~~~~~~

- Robert T. McGibbon
- Martin K. Scherer
- Alex Yang +
- Fabian Paul
- Kirill Shmilovich +
- Lucian Krapp +
- Sander Roet +
- David W.H. Swenson
- Ray A. Matsumoto
- Jack Greisman
- Marius van Niekerk +
- Patrick Kunzmann +
- Matthew W. Thompson
- Justin R. Porter
- Richard Banh +
- sefalkner +

A total of 16 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.


v1.9.3 (May 17, 2019)
---------------------

- Fix CI (#1416, #1420, #1444)
- Add compute_center_of_geometry (#1405)
- Fix a test failure in test_reporters.py caused by merge of #1431 (#1443)

- Reporters no longer override user request to write unit cell information (#1431)
- Add XTCReporter for OpenMM (#1403)
- [xtc] Fix bugs in striding with atom_indices and seek+stride (#1449)

- Avoid infinite recursion error in mol2 parser (#1426)
- [formats/mol2] add more checks to element parsing (#1407)
- Replace strip() with split() in `mol2.py` (#1378)

- Use and set resSeq attribute in Topology.to_openmm() and from_openmm() (#1424)
- fix parallel reduction error (#1419)
- Fixes 'Buffer dtype mismatch' error on 64-bit Windows (#1409)

- add RMSF analysis (#1414)
- allow RMSD calls when ref_atom_indices and atom_indices are used (#1392)
- Notebook tests: `from __future__` must come first (#1401)

- [setup] do not enforce clang/std++ on osx (#1400)
- silence cython related numpy warnings (#1391)
- Prep py37, some bugfixes (#1388)
- Ensure 'bond_value' is a string (#1382)
- fix typo in docs (#1381)


Authors
~~~~~~~

- Carlos Hernández
- John Chodera
- Jack Greisman
- jgilaber
- Sunhwan Jo
- Ray A. Matsumoto
- Robert T. McGibbon
- João Rodrigues
- Shyam Saladi
- Martin K. Scherer
- David W.H. Swenson
- Matthew W. Thompson
- Lee-Ping Wang

A total of 12 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.


v1.9.2 (July 30, 2018)
----------------------
We're please to announce the release of MDTraj 1.9.2. This version has a number of bug fixes and improvements for trajectory parsing and conversion.


- Fix bug in TINKER ARC reader (#1371)
- Improved mdconvert error message (#1368)
- Striding relative to current position in XTC and TRR (#1364)
- Return last successful read frame for DCD (#1358)
- Handle stride like numpy for DCDs (#1352)
- Fix pickling of virtual site's element field (#1350)
- Compile geometry extension with OpenMP (#1349)
- Ensure correct dtype in neighborlist box vectors (#1344)
- Added support for prm7 topology file extension (#1334)
- Added efficient stride handling fo TRR (#1332)
- Use byte offsets between frames for stride of XTCs (#1331)
- Updated the calculation of chi5 (#1322, #1323)
- Added testing against conda-forge channel (#1310)
- Port [OpenMM bond order](https://github.com/pandegroup/openmm/pull/1668) representation into MDTraj. Implements the `Bond` class to Topology and updates the Mol2 reader to use bond_order field (#1308)

Authors
~~~~~~~

- Carlos Hernández
- Guillermo Pérez-Hernández
- Matthew Harrigan
- Lester Hedges +
- Robert T. McGibbon
- Levi Naden +
- Fabian Paul
- Justin R. Porter
- Martin K. Scherer
- Xianqiang Sun +
- David W.H. Swenson +
- Lee-Ping Wang

A total of 11 people contributed to this release.
People with a "+" by their names contributed a patch for the first time.


v1.9 (September 3, 2017)
------------------------

- [xtc] ``approx_nframes`` returns at least one (#1265)
- Make ``compute_directors`` user-facing (#1260)
- Add differentiable contacts option (#1247)
- Remove link to forum (#1237)
- Skip renumbering if no bonds in mol2 (#1238)
- Add a bunch of Van Der Waals values (#1174)
- [geometry] Fix compatibility with old visual studio for Python 2.7 (#1233)
- Implement ``compute_average_structure`` (#1221)
- Fix import of ``load_stk`` (#1231)
- Fix bugs in load with atom_indices and frame args (#1227)
- Fix conda test running (#1228)
- Amber restart file convention (#1223)
- Install path for zlib on linux too (#1208)
- Fix transform calculation and Transform object to be more general (#1254)
- Add O1 as alternative definition for water oxygen (#1257)
- Fix precentering overflow error in center.c (#1283)
- Add chi5 angle computation (#1291)
- Fix the build bug caused by incorrect plumbing of the numpy include path
- into ``cflags`` (#1290)
- Make RDF ``pairs`` argument required (#1288)
- Refresh tests (#1266)
- Remove PyPI downloads badge (#1293)
- Extracting velocities/forces from TRR files (hidden API) (#1294)
- Add "in" selection to selection language (#1268)
- Handle a single frame being passed to sparta+ (#1295)

v1.8 (November 9, 2016)
-----------------------

- PR #1202: ``mdtraj.html`` has been removed. We recommend using
  ``nglview`` for visualizing MDTraj trajectory objects.
- PR #1204: Fix search functionality with docs
- PR #1167: Fix corner case in distancekernel.h
- PR #1190: Fix issue with rmsd precentered = True and atom_indices != None
- PR #1106: Speed up image_molecules
- PR #1182: Add 'sidechain' and 'sidechain-heavy' options to compute_contacts
- PR #1180: Handle unexpected keyword arguments gracefully in psf and prmtop parsers
- PR #1171: Remove unnecessary restriction on iterload
- PR #1170: Load single-element path lists without a copy
- PR #1165: There should never be zero bins in Voxels class
- PR #1158: Update deprecated use of scipy.stats.nanmean
- PR #1153: [formats/XTC] in case of an out of bounds seek, raise IOError
- PR #1161: Fix typos in examples
- PR #1130: Automatically test examples to make sure they work
- PR #1155: Update wording for simulation-with-openmm.ipynb
- PR #1146: Ensure box vectors have right dtype
- PR #1145: Check that file exists before trying to open it
- PR #1139: Optimize baker_hubbard and wernet_nilsson functions
- PR #1137: Allow standard_names as a keyword argument to md.load()
- PR #1132: Fix bug in hoomdxml reader
- PR #1125: Support Gromacs TNG files
- PR #1123: Add md.join(trajs)

v1.7.2 (May 2, 2016)
--------------------

- Small fix to developer tools so docs get uploaded.

v1.7 (May 2, 2016)
------------------

We're please to announce the release of MDTraj 1.7. In addition to the
usual fixes and improvements, MDTraj has gained the ability to image
molecules in trajectories. So far, it's worked very well even on
complicated systems like multi-molecule proteins. Look forward to future
enhancements to this new feature! Some other highlights include:

- New ``compute_neighborlist()`` function (#1057)
- Add option to skip standardization of atom and residue names during
  ``load_pdb`` (#1061)
- Function for imaging molecules (#1058)
- New optional argument ``periodic`` for ``compute_contacts`` (#1072)
- Refresh documentation (#1067, #1074, #1075)
- Rewrite geometry code in modern c++ (#1077)
- Fix issue with ``Topoplogy.from_openmm`` (#1089)


v1.6 (February 15, 2016)
------------------------

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

.. vim: tw=75
