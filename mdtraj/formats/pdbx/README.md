pdbx
====

OpenMM PDBx/mmCIF Reader adapted for MDTraj, with the following changes:

* **MDTraj-native Topology methods**

  * Replaced all `Topology.addChain/addResidue/addAtom` calls and generator methods (`chains()`, `residues()`, `atoms()`) with MDTraj’s properties and methods (`add_chain`, `add_residue`, `add_atom` and `.chains`, `.residues`, `.atoms`).

* **NumPy-based positions (no OpenMM units/Vec3)**

  * Removed OpenMM’s `Vec3` and `Quantity` dependencies. All coordinates are stored as raw NumPy arrays in nanometers. Positions stored as numpy arrays by default (removed `get_positions(as_numpy)`).


* **Bond construction & deduplication**

  * After calling `topology.create_standard_bonds()`, the code reads both `struct_conn` and `chem_comp_bond` categories. It then removes duplicate bonds, merging type/order with warnings if they conflict.

* **Atom/residue name‐replacement tables**

  * Adopted atom and residue name replacements consistent with MDTraj's existing PDB reader.

* **Box‐lengths/angles handling**

  * Instead of using OpenMM’s `computeLengthsAndAngles`, we manually read `_cell.length_{a,b,c}` (in Ångström → convert to nm) and `_cell.angle_{α,β,γ}` (in degrees).
  * Moved periodic box vectors and angles handling from `Topology` into `Trajectory`.
  * It is consistent with mdtraj `Trajectory` expectations.

* **B‐factor support in `pdbxfile`**

  * `writeModel` now accepts an optional `bfactors` list (defaults to 0.0 if `None`).

* **MDTraj version stamping**

  * The header prints `# Created with mdtraj <version>, <date>`.

* **Docstrings in MDTraj style**

  * Every public method now has a NumPy‐style docstring (parameters, returns, exceptions).

* **Exception specificity**

  * Wherever the original code used `except:` or caught all exceptions, we now catch only `IndexError`, `KeyError`, `ValueError`, or `AttributeError`.

* **Removal of `__future__` imports**

* **Replaced `xrange` with `range`**

* **Converted `%`‐style formatting to f-strings**

* **Expanded imports onto separate lines**

* **Fixes parsing cif resid as integer in `pdbxfile`**

* **Linting**

