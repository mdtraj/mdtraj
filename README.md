## Molecular Dynamics Trajectory: File IO and manipulation

[![Build Status](https://travis-ci.org/rmcgibbo/mdtraj.png)](https://travis-ci.org/rmcgibbo/mdtraj)

This library supports the reading and writing of molecular dynamics trajectories in a variety of formats. Currently, there is full support for [PDB](http://deposit.rcsb.org/adit/docs/pdb_atom_format.html), [DCD](https://www-s.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html),
[XTC](http://manual.gromacs.org/online/xtc.html), [TRR](http://www.gromacs.org/Documentation/File_Formats/.trr_File),
[binpos](https://www-s.ks.uiuc.edu/Research/vmd/plugins/molfile/binposplugin.html), [AMBER NetCDF](http://ambermd.org/netcdf/nctraj.html),
[AMBER mdcrd](http://ambermd.org/formats.html), and [MDTraj HDF5](https://github.com/rmcgibbo/mdtraj/wiki/HDF5-Trajectory-Format).
There is partial support for the TINKER arc format (only reading).

MDTraj is structured around a convent numpy-based Trajectory object. Its RMSD (optimal cartesian
root-mean-square deviation) library gives you access to fast parallel SSE3-based structural
deviations using the quaternion characteristic polynomial (Theobald QCP) method at 4x the
speed of the original Theobald code and over 3x as fast as Theobald code modified to use GotoBLAS.

For details, see the [documentation](http://rmcgibbo.github.io/mdtraj/).
