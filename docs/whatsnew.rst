**********
What's New
**********

These are new features and improvements of note in each release.

v0.7.0 (February 21, 2014)
--------------------------
- New geometry functions `md.compute_contacts` and `md.geometry.squareform` for residue-residue contact maps
- Fix segfault in `md.rmsd` when using the optional `atom_indices` kwarg
- `md.compute_phi`, `md.compute_psi`, and `md.compute_omega` now return the correct atom indices, as their docstring always said.
- Topology `Element` instances are now properly immutable
- Small bugfixes to `baker_hubbard`, and better docstring
- Automatic installation of `pandas` and `simtk.unit` via setuptools' `install_requires`.
- Small bugfix to mdcrd loading with stride
- `superpose` now correctly translates the final structure, and doesn't recenter the reference structure

v0.6.1 (February 11, 2014)
--------------------------
- `Trajectory.join(discard_overlapping_frames=True)` is criterion for detecting overlappign frames is more realistic
- We now support installation via conda, and are supplying conda binaries
- `md.load()` is much faster when loading multiple trajectory files
- Bugfixes for pandas 0.13.0 release, detection of zinc atoms in PDB files
- Geometry functions are more resiliant to segfaults from bad user parameters
- Fix intermittent RMSD segfaults from invalid memory access
- Fix RMSD centering bug with memory alignment after restrict_atoms

v0.6.0 (January 21, 2014)
-------------------------
- `md.rmsd()` signature changed to be more understandable
- All file objects now have a `__len__` function.
- Small bugfixes related to vsites.

v0.5.1 (January 4, 2014)
------------------------
- Minor bug fix when no dihedrals match specification
- Add __str__ to Topology parts
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
- Hydrogen bond identification with baker_hubbard
- Rotation/translation to superpose trajectories
- New RMSD API. It's much simpler and much more memory efficient
- Full support for computing all of the chi angles
- Add seek/tell methods to all of the trajectory file objects
- New top level memory efficient iterload method for chunked trajectory loading
