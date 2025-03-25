pdbx
====

OpenMM PDBx/mmCIF Reader adapted for MDTraj, with the following changes:

- Removed dependency on `openmm.units`; positions are now numpy arrays without units.
- Removed `Vec3`; all vector operations use numpy arrays directly.
- Positions stored as numpy arrays by default (removed `get_positions(as_numpy)`).
- Updated documentation format to align with MDTraj conventions.
- Replaced OpenMM's `Topology` class with MDTraj's native `Topology`, converting methods (e.g., `chains()`, `residues()`, `atoms()`) into properties (`chains`, `residues`, `atoms`) to align with MDTraj's API style.
- Moved periodic box vectors and angles handling from `Topology` into `Trajectory`.
- Adopted atom and residue name replacements consistent with MDTraj's existing PDB reader for uniform trajectory handling.