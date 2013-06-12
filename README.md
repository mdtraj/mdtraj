## Molecular Dynamics Trajectory: File IO and manipulation

[![Build Status](https://travis-ci.org/rmcgibbo/mdtraj.png)](https://travis-ci.org/rmcgibbo/mdtraj)

This library supports the reading and writing of molecular dynamics trajectories in a variety of formats.

Currently, there is full support for reading and writing from
 - [PDB](http://deposit.rcsb.org/adit/docs/pdb_atom_format.html)
 - [DCD](https://www-s.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html)
 - [XTC](http://manual.gromacs.org/online/xtc.html)
 - [TRR](http://www.gromacs.org/Documentation/File_Formats/.trr_File)
 - [binpos](https://www-s.ks.uiuc.edu/Research/vmd/plugins/molfile/binposplugin.html)
 - [AMBER NetCDF](http://ambermd.org/netcdf/nctraj.html)
 - [MDTraj HDF5](https://github.com/rmcgibbo/mdtraj/wiki/HDF5-Trajectory-Format)

The central data structure is the `Trajectory` object, which contains the positions from mutliple frames
of a molecular dynamics simulations in a numpy array and the system topology. It also can contain
information about the unitcell geometry. The trajectory object is designed to be convenient to work with.
