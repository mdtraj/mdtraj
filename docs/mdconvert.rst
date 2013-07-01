Command-line trajectory conversion: ``mdconvert``
=================================================

``mdconvert`` is a command-line script installed with MDTraj to convert
molecular dynamics trajectories between formats. The DCD, XTC, TRR,
binpos, NetCDF, binpos, and HDF5 formats are supported (.xtc, .nc, .trr, .h5,
.pdb, .binpos, .dcd). ``mdconvert`` is memory-efficient, and processes
trajectories in a chunked, streaming fashion. It is capable of converting
trajectory files which cannot be fully loaded into memory. It can also
concatenate trajectories, convert only a subset of the atoms in a trajectory
(i.e. strip solvent molecules), and down-sample trajectories by extract only a
subset of the frames.

After installing the library, it should be in your  ``$PATH``. You can check
this from the command line with this command. ::
  
  $ which mdconvert
  

Here's the ``mdconvert`` help text. ::

  $ mdconvert -h
  usage: mdconvert [-h] -o OUTPUT [-c CHUNK] [-f] [-s STRIDE] [-i INDEX]
                   [-a ATOM_INDICES] [-t TOPOLOGY]
                   input [input ...]

  Convert molecular dynamics trajectories between formats. The DCD, XTC, TRR,
  binpos, NetCDF, binpos, and HDF5 formats are supported (.xtc, .nc, .trr, .h5,
  .pdb, .binpos, .dcd)

  positional arguments:
    input                 path to one or more trajectory files. Multiple
                          trajectories, if supplied, will be concatenated
                          together in the output file in the order supplied. all
                          of the trajectories should be in the same format. the
                          format will be detected based on the file extension

  required arguments:
    -o OUTPUT, --output OUTPUT
                          path to the save the output. the output format will
                          chosen based on the file extension (.xtc, .nc, .trr,
                          .h5, .pdb, .binpos, .dcd)

  optional arguments:
    -h, --help            show this help message and exit
    -c CHUNK, --chunk CHUNK
                          number of frames to read in at once. this determines
                          the memory requirements of this code. default=1000
    -f, --force           force overwrite if output already exsits
    -s STRIDE, --stride STRIDE
                          load only every stride-th frame from the input
                          file(s), to subsample.
    -i INDEX, --index INDEX
                          load a *specific* set of frames. flexible, but
                          inefficient for a large trajectory. specify your
                          selection using (pythonic) "slice notation" e.g. '-i
                          N' to load the the Nth frame, '-i -1' will load the
                          last frame, '-i N:M to load frames N to M, etc. see
                          http://bit.ly/143kloq for details on the notation
    -a ATOM_INDICES, --atom_indices ATOM_INDICES
                          load only specific atoms from the input file(s).
                          provide a path to file containing a space, tab or
                          newline separated list of the (zero-based) integer
                          indices corresponding to the atoms you wish to keep.
    -t TOPOLOGY, --topology TOPOLOGY
                          path to a PDB file. this will be used to parse the
                          topology of the system. it's optional, but useful. if
                          specified, it enables you to output the coordinates of
                          your dcd/xtc/trr/netcdf/binpos as a PDB file. If
                          you're converting *to* .h5, the topology will be
                          stored inside the h5 file.