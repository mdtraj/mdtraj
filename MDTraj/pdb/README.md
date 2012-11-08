openmmpdb
=========

OpenMM PDB Reader and Topology with a few changes

- Removed all of the units, so that this can operate without simtk.units
- Removed Vec3, and just used numpy arrays
- Made positions be numpy by default (no get_positions(as_numpy))
- Added to_bytearray and from_bytearray to topology, so that it can serialize itself
