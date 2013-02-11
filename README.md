## Molecular Dynamics Trajectory: File IO and manipulation

[![Build Status](https://travis-ci.org/rmcgibbo/mdtraj.png)](https://travis-ci.org/rmcgibbo/mdtraj)

Reads from and writes  PDB, DCD, XTC, BinPos, and PyTables hdf5 files. This repository is a working sketch of a replacement for MSMBuilder.Trajectory, and may be directly used in https://github.com/tjlane/odin

The Trajectory data structure contains a numpy array (n_frames x n_atoms x 3), an OpenMM style object oriented topology
(Residue, Chain, Atom, etc). It also contains the timestep, and the topology should also store the box dimensions.
