## Molecular Dynamics Trajectory: File IO and manipulation

[![Build Status](https://travis-ci.org/rmcgibbo/mdtraj.png)](https://travis-ci.org/rmcgibbo/mdtraj)

Reads from and writes  PDB, DCD, XTC, BinPos, and PyTables hdf5 files. This repository is a working sketch of a replacement for MSMBuilder.Trajectory, and may be directly used in https://github.com/tjlane/odin

The Trajectory data structure contains a numpy array (n_frames x n_atoms x 3), an OpenMM style object oriented topology
(Residue, Chain, Atom, etc). It also contains the timestep, and the topology should also store the box dimensions.

## HDF5 Compression
Within pytables, we use by default a lossy 16bit fixed-width encoding for the atom positions, in nm. The mean
absolute error of the compression is roughly ~5e-4 nm, which is acceptible. Only numbers between -32 and +32
nanometers can be encoded with this scheme. The performance, as compared to a standard IEEE float16 compression,
is shown in this IPython notebook: http://nbviewer.ipython.org/5399019
