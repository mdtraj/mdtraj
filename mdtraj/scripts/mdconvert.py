##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
# 
# Authors: Robert McGibbon
# Contributors:
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################


###############################################################################
# Imports
###############################################################################
from __future__ import print_function

import os
import sys
import glob
import warnings
import functools
import operator
from argparse import ArgumentParser

import numpy as np
import mdtraj as md
from mdtraj.core.trajectory import _parse_topology
from mdtraj.utils import in_units_of
from mdtraj.utils.six import iteritems

###############################################################################
# Crappy class that should go elsewhere
###############################################################################

###############################################################################
# Globals
###############################################################################

formats = {'.dcd': md.formats.DCDTrajectoryFile,
           '.xtc': md.formats.XTCTrajectoryFile,
           '.trr': md.formats.TRRTrajectoryFile,
           '.binpos': md.formats.BINPOSTrajectoryFile,
           '.nc': md.formats.NetCDFTrajectoryFile,
           '.netcdf': md.formats.NetCDFTrajectoryFile,
           '.h5': md.formats.HDF5TrajectoryFile,
           '.lh5': md.formats.LH5TrajectoryFile,
           '.pdb': md.formats.PDBTrajectoryFile}

fields = {'.trr': ('xyz', 'time', 'step', 'box', 'lambda'),
          '.xtc': ('xyz', 'time', 'step', 'box'),
          '.dcd': ('xyz', 'cell_lengths', 'cell_angles'),
          '.nc': ('xyz', 'time', 'cell_lengths', 'cell_angles'),
          '.netcdf': ('xyz', 'time', 'cell_lengths', 'cell_angles'),
          '.binpos': ('xyz',),
          '.lh5': ('xyz', 'topology'),
          '.h5': ('xyz', 'time', 'cell_lengths', 'cell_angles',
                  'velocities', 'kineticEnergy', 'potentialEnergy',
                  'temperature', 'lambda', 'topology'),
          '.pdb': ('xyz', 'topology', 'cell_angles', 'cell_lengths')}

units = {'.xtc': 'nanometers',
         '.trr': 'nanometers',
         '.binpos': 'angstroms',
         '.nc': 'angstroms',
         '.netcdf': 'angstroms',
         '.dcd': 'angstroms',
         '.h5': 'nanometers',
         '.lh5': 'nanometers',
         '.pdb': 'angstroms'}

###############################################################################
# Utility Functions
###############################################################################


ext = lambda fn: os.path.splitext(fn)[1]


class _Warner(object):
    def __init__(self):
        self.active = True

    def __call__(self, msg):
        if self.active:
            print('Warning:', msg, file=sys.stderr)
warn = _Warner()


def index(str):
    if str.count(':') == 0:
        return int(str)

    elif str.count(':') == 1:
        start, end = [(None if e == '' else int(e)) for e in str.split(':')]
        step = None
    elif str.count(':') == 2:
        start, end, step = [(None if e == '' else int(e)) for e in str.split(':')]

    return slice(start, end, step)


###############################################################################
# Code
###############################################################################


def parse_args():
    """Parse the command line arguments and perform some validation on the
    arguments

    Returns
    -------
    args : argparse.Namespace
        The namespace containing the arguments
    """
    extensions = ', '.join(list(formats.keys()))
    parser = ArgumentParser(description='''Convert molecular dynamics
    trajectories between formats. The DCD, XTC, TRR, PDB, binpos, NetCDF,
    binpos, LH5, and HDF5 formats are supported (%s)''' % extensions)
    parser.add_argument('input', nargs='+', help='''path to one or more
                    trajectory files. Multiple trajectories, if supplied, will
                    be concatenated together in the output file in the order
                    supplied. all of the trajectories should be in the same
                    format. the format will be detected based on the file
                    extension''')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-o', '--output', required=True,
                          help='''path to the save the output. the output
                          format will chosen based on the file extension
                          (%s)''' % extensions)
    # dirty hack to move the 'optional arguments' group to the end. such that
    # the 'required arguments' group shows up before it.
    parser._action_groups.append(parser._action_groups.pop(1))
    parser.add_argument('-c', '--chunk', default=1000, type=int,
                        help='''number of frames to read in at once. this
                        determines the memory requirements of this code.
                        default=1000''')
    parser.add_argument('-f', '--force', action='store_true',
                        help='''force overwrite if output already exsits''')
    parser.add_argument('-s', '--stride', default=1, type=int, help='''load
                        only every stride-th frame from the input file(s),
                        to subsample.''')
    parser.add_argument('-i', '--index', type=index, help='''load a *specific*
                        set of frames. flexible, but inefficient for a large
                        trajectory. specify your selection using (pythonic)
                        "slice notation" e.g. '-i N' to load the the Nth
                        frame, '-i -1' will load the last frame, '-i N:M to
                        load frames N to M, etc. see http://bit.ly/143kloq
                        for details on the notation''')
    parser.add_argument('-a', '--atom_indices',  type=str,
                        help='''load only specific atoms from the input file(s).
                        provide a path to file containing a space, tab or
                        newline separated list of the (zero-based) integer
                        indices corresponding to the atoms you wish to keep.''')
    parser.add_argument('-t', '--topology', type=str, help='''path to a
                        PDB/prmtop file. this will be used to parse the topology
                        of the system. it's optional, but useful. if specified,
                        it enables you to output the coordinates of your
                        dcd/xtc/trr/netcdf/binpos as a PDB file. If you\'re
                        converting *to* .h5, the topology will be stored
                        inside the h5 file.''')

    args = parser.parse_args()

    if not args.force and os.path.exists(args.output):
        parser.error('file exists: %s' % args.output)

    # rebuild the input list, doing any glob expansions
    # necessary
    input = []
    for fn in args.input:
        if not os.path.exists(fn):
            if '*' in fn:
                input.extend(glob.glob(fn))
            else:
                parser.error('No such file: %s' % fn)
        elif os.path.isdir(fn):
            parser.error('%s: Is a directory' % fn)
        elif not os.path.isfile(fn):
            parser.error('%s: Is not a file' % fn)
        else:
            input.append(fn)
    args.input = input

    for fn in args.input:
        if not ext(fn) in formats:
            parser.error("%s: '%s' is not a known extension" % (fn, ext(fn)))

    extensions = list(map(ext, args.input))
    if any(e != extensions[0] for e in extensions):
        parser.error("all input trajectories do not have the same extension")

    if not ext(args.output) in formats:
        parser.error("%s: '%s' is not a known extension" % (args.output,
                     ext(args.output)))

    if args.atom_indices is not None and not os.path.isfile(args.atom_indices):
        parser.error('no such file: %s' % args.atom_indices)

    if args.stride <= 0:
        parser.error('stride must be positive')
    if args.chunk <= 0:
        parser.error('chunk must be positive')

    if args.index and len(args.input) > 1:
        parser.error('index notation only allowed with a single input trajectory')
    if args.index and args.stride != 1:
        parser.error('stride and index selections are incompatible')
    if args.index is not None:
        args.chunk = None

    if args.topology is not None and not os.path.isfile(args.topology):
        parser.error('no such file: %s' % args.topology)

    if ((args.topology is None and not all(ext(e) in ['.h5', '.lh5', '.pdb'] for e in args.input))
                and ext(args.output) in ['.h5', '.lh5', '.pdb']):
        parser.error('to output a %s file, you need to supply a topology (-t, or --topology)' % ext(args.output))

    if args.chunk is not None and (args.chunk % args.stride != 0):
        parser.error('--stride must be a divisor of --chunk')

    return args


def main(args, verbose=True):
    """Run the main script.

    Parameters
    ----------
    args : argparse.Namespace
        The collected command line arguments
    """
    if args.atom_indices is not None:
        atom_indices = np.loadtxt(args.atom_indices, int)
    else:
        atom_indices = None

    out_x = ext(args.output)
    out_units = units[out_x]
    out_fields = fields[out_x]
    OutFileFormat = formats[out_x]

    in_x = ext(args.input[0])
    InFileFormat = formats[in_x]

    if args.topology is not None:
        topology = _parse_topology(args.topology)
    else:
        topology = None

    if topology is not None and atom_indices is not None:
        topology = topology.subset(atom_indices)

    n_total = 0
    if args.index is not None:
        assert len(args.input) == 1
        # when chunk is None, we load up ALL of the frames. this isn't
        # strictly necessary, and it costs more memory, but it's ALOT
        # harder to get the code correct when we need to use data[start:end]
        # notation when all of the data isn't loaded up at once. it's easy
        # for hdf5 and netcdf, but for the others...
        assert args.chunk is None

    # this is the normal invocation pattern, but for PDBTrajectoryFile it's
    # different
    outfile_factory = functools.partial(OutFileFormat, args.output, 'w',
                        force_overwrite=args.force)

    with outfile_factory() as outfile:
        for fn in args.input:
            assert in_x == ext(fn)
            with InFileFormat(fn, 'r') as infile:

                while True:
                    data, in_units, n_frames = read(infile, args.chunk, stride=args.stride,
                                                    atom_indices=atom_indices)
                    if n_frames == 0:
                        break

                    if topology is not None:
                        # if the user supplied a topology, we should probably
                        # do some simple checks
                        if data['xyz'].shape[1] != topology._numAtoms:
                            warnings.warn('sdsfsd!!!!')
                        data['topology'] = topology

                    # if they want a specific set of frames, get those
                    # with slice notation
                    if args.index is not None:
                        _data = {}
                        for k, v in iteritems(data):
                            if isinstance(v, np.ndarray):
                                # we don't want the dimensionality to go deficient
                                if isinstance(args.index, int):
                                    _data[k] = v[np.newaxis, args.index]
                                else:
                                    _data[k] = v[args.index]
                            elif isinstance(v, md.Topology):
                                _data[k] = v
                            else:
                                raise RuntineError()
                        data = _data
                        print(list(data.keys()))
                        n_frames = len(data['xyz'])

                    convert(data, in_units, out_units, out_fields)
                    write(outfile, data)
                    n_total += n_frames

                    if verbose:
                        sys.stdout.write('\rconverted %d frames, %d atoms' % (n_total, data['xyz'].shape[1]))
                        sys.stdout.flush()

    if verbose:
        print(' ')


def write(outfile, data):
    """Write data out to a file

    This is a small wrapper around the native write() method on the
    XXXTRajectoryFile objects that is necessary to make sure we pass the
    right arguments in the right position

    Parameters
    ----------
    outfile : TrajectoryFile
        An open trajectory file with a write() method
    data : dict
        A dict with the data to write in it.
    """
    if isinstance(outfile, md.formats.XTCTrajectoryFile):
        outfile.write(data.get('xyz', None), data.get('time', None),
                      data.get('step', None), data.get('box', None))

    elif isinstance(outfile, md.formats.TRRTrajectoryFile):
        outfile.write(data.get('xyz', None), data.get('time', None),
                      data.get('step', None), data.get('box', None),
                      data.get('lambd', None))

    elif isinstance(outfile, md.formats.DCDTrajectoryFile):
        outfile.write(data.get('xyz', None), data.get('cell_lengths', None),
                      data.get('cell_angles', None))

    elif isinstance(outfile, md.formats.BINPOSTrajectoryFile):
        outfile.write(data.get('xyz', None))

    elif isinstance(outfile, md.formats.PDBTrajectoryFile):
        lengths, angles = None, None
        for i, frame in enumerate(data.get('xyz')):
            if 'cell_lengths' in data:
                lengths = data['cell_lengths'][i]
            if 'cell_angles' in data:
                angles = data['cell_angles'][i]
                
            outfile.write(frame, data.get('topology', None), i, lengths, angles)

    elif isinstance(outfile, md.formats.NetCDFTrajectoryFile):
        outfile.write(data.get('xyz', None), data.get('time', None),
                      data.get('cell_lengths', None), data.get('cell_angles', None))

    elif isinstance(outfile, md.formats.HDF5TrajectoryFile):
        outfile.write(data.get('xyz', None), data.get('time', None),
                      data.get('cell_lengths', None), data.get('cell_angles', None),
                      data.get('velocities', None), data.get('kineticEnergy', None),
                      data.get('potentialEnergy', None), data.get('temperature', None),
                      data.get('lambda', None))
        if outfile.topology is None:
            # only want to write the topology once if we're chunking
            outfile.topology = data.get('topology', None)

    elif isinstance(outfile, md.formats.LH5TrajectoryFile):
        outfile.write(data.get('xyz', None))
        if outfile.topology is None:
            # only want to write the topology once if we're chunking
            outfile.topology = data.get('topology', None)
    else:
        raise RuntimeError()


def read(infile, chunk, stride, atom_indices):
    """Read data from the infile.

    This is a small wrapper around the read() method on the XXXTrajectoryFile
    that performs the read and then puts the results in a little dict. It also
    returns the distance units that the file uses.
    """
    if not isinstance(infile, md.formats.PDBTrajectoryFile):
        _data = infile.read(chunk, stride=stride, atom_indices=atom_indices)

    if isinstance(infile, md.formats.PDBTrajectoryFile):
        if infile.closed:
            # signal that we're done reading this pdb
            return None, None, 0

        if atom_indices is None:
            atom_indices = slice(None)
            topology = infile.topology
        else:
            topology = infile.topology.subset(atom_indices)

        data = {'xyz': infile.positions[::stride, atom_indices, :],
                'topology': topology}
        if infile.unitcell_lengths is not None:
            data['cell_lengths'] =np.array([infile.unitcell_lengths] * len(data['xyz']))
            data['cell_angles'] = np.array([infile.unitcell_angles] * len(data['xyz']))
        in_units = 'angstroms'
        infile.close()

    elif isinstance(infile, md.formats.XTCTrajectoryFile):
        data = dict(zip(fields['.xtc'], _data))
        in_units = 'nanometers'

    elif isinstance(infile, md.formats.TRRTrajectoryFile):
        data = dict(zip(fields['.trr'], _data))
        in_units = 'nanometers'

    elif isinstance(infile, md.formats.DCDTrajectoryFile):
        data = dict(zip(fields['.dcd'], _data))
        in_units = 'angstroms'

    elif isinstance(infile, md.formats.BINPOSTrajectoryFile):
        data = {'xyz': _data}
        in_units = 'angstroms'

    elif isinstance(infile, md.formats.NetCDFTrajectoryFile):
        data = dict(zip(fields['.nc'], _data))
        in_units = 'angstroms'

    elif isinstance(infile, md.formats.HDF5TrajectoryFile):
        data = dict(zip(fields['.h5'], _data))
        data['topology'] = infile.topology  # need to hack this one in manually
        if atom_indices is not None:
            data['topology'] = data['topology'].subset(atom_indices)
        in_units = 'nanometers'

    elif isinstance(infile, md.formats.LH5TrajectoryFile):
        data = {'xyz': _data}
        data['topology'] = infile.topology  # need to hack this one in manually
        if atom_indices is not None:
            data['topology'] = data['topology'].subset(atom_indices)
        in_units = 'nanometers'

    else:
        raise RuntimeError

    data = dict((k, v) for k, v in data.items() if v is not None)
    return data, in_units, (0 if 'xyz' not in data else len(data['xyz']))


def convert(data, in_units, out_units, out_fields):
    # do unit conversion
    if 'xyz' in out_fields and 'xyz' in data:
        data['xyz'] = in_units_of(data['xyz'], in_units, out_units, inplace=True)
    if 'box' in out_fields:
        if 'box' in data:
            data['box'] = in_units_of(data['box'], in_units, out_units, inplace=True)
        elif 'cell_angles' in data and 'cell_lengths' in data:
            a, b, c = data['cell_lengths'].T
            alpha, beta, gamma = data['cell_angles'].T
            data['box'] = np.dstack(md.utils.unitcell.lengths_and_angles_to_box_vectors(a, b, c, alpha, beta, gamma))
            data['box'] = in_units_of(data['box'], in_units, out_units, inplace=True)
            del data['cell_lengths']
            del data['cell_angles']

    if 'cell_lengths' in out_fields:
        if 'cell_lengths' in data:
            data['cell_lengths'] = in_units_of(data['cell_lengths'], in_units, out_units, inplace=True)
        elif 'box' in data:
            a, b, c, alpha, beta, gamma = md.utils.unitcell.box_vectors_to_lengths_and_angles(data['box'][:, 0], data['box'][:, 1], data['box'][:, 2])
            data['cell_lengths'] = np.vstack((a, b, c)).T
            data['cell_angles'] = np.vstack((alpha, beta, gamma)).T
            data['cell_lengths'] = in_units_of(data['cell_lengths'], in_units, out_units, inplace=True)
            del data['box']

    ignored_keys = ["'%s'" % s for s in set(data) - set(out_fields)]
    formated_fields = ', '.join("'%s'" % o for o in out_fields)
    if len(ignored_keys) > 0:
        warn('%s data from input file(s) will be discarded. '
             'output format only supports fields: %s' % (', '.join(ignored_keys),
                                                         formated_fields))
        warn.active = False

    return data

def entry_point():
    args = parse_args()
    main(args)
    
if __name__ == '__main__':
    entry_point()
