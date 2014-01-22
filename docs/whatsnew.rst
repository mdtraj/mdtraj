**********
What's New
**********

These are new features and improvements of note in each release.

v0.6.0 (January 21, 2014)
-----------------
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
