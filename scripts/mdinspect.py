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

"""Sanity checking for molecular dynamics trajectories. This script is currently
a work in progress. Contributions are encouraged. 
"""
#------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------

from __future__ import print_function
import os
import sys
import warnings
import functools
import operator
from argparse import ArgumentParser

import numpy as np
import mdtraj as md
from mdtraj.utils import import_, ilen
from mdtraj.geometry.internal import COVALENT_RADII

spatial = import_('scipy.spatial')

#------------------------------------------------------------------------------
# Code
#------------------------------------------------------------------------------

class NoTopologyError(Exception):
    def __init__(self):
        super(NoTopologyError, self).__init__("One more more of the "
            "trajectory files should contain topology information (i.e. "
            "either HDF5 or PDB)")


def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='+', help='''Input trajectory file(s),
        in any supported format. The files are assumed to be different
        trajectories for the same system. One or more of the trajectory
        files should contain topology information (i.e. one of the files
        should be either a PDB or HDF5)''')
    #parser.add_argument('-n', '--noload', action='store_true', help='''Do not load the coordinate data from the trajectory, for example if the trajectory is too large to load into memory. Only a limited number of checks will be done.''')
    parser.add_argument('--bond-low', type=float, help='''Minimum fraction of sum of covalent radii for bonded atoms. Default=0.4''', default=0.4)
    parser.add_argument('--bond-high', type=float, help='''Maximum fraction of sum of covalent radii for bonded atoms. Default=1.2''', default=1.2)
    return parser.parse_args(), parser


def main(args, parser):
    inspector = Inspector(args.bond_low, args.bond_high)
    for f in sorted(args.files, cmp=cmp_pdb):
        if not os.path.exists(f):
            parser.error("File '%s' does not exist" % f)
        if not os.path.isfile(f):
            parser.error("File '%s' is not a file" % f)

        try:
            inspector.load(f)
        except NoTopologyError as e:
            parser.error(str(e))


def cmp_pdb(a, b):
    """String comparision function, for sorting, that puts things
    ending in .pdb at the beginning"""
    if a.endswith('.pdb') and b.endswith('.pdb'):
        return 0
    if a.endswith('.pdb'):
        return -1
    if b.endswith('.pdb'):
        return 1
    if a < b:
        return -1
    if a > b:
        return 1
    return 0


class Inspector(object):
    def __init__(self, bond_low, bond_high):
        self._printed_section = False
        self.t = None
        self.fn = None
        self.top = None

        self.bond_low = bond_low
        self.bond_high = bond_high

    def load(self, fn):
        self.fn = fn
        ext = os.path.splitext(fn)[1]

        if self.top is not None or ext in ['.pdb', '.h5']:
            self.load_topology(fn)
        else:
            self.load_no_topology(fn)

    def load_topology(self, fn):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.t = md.load(fn, top=self.top)
            self.top = self.t.topology

        self.check_unitcell()
        self.check_topology()
        self.check_bonds()
        self.check_positions()

    def load_no_topology(self, fn):
        raise NoTopologyError()

    def check_unitcell(self):
        self.section('unitcell')
        if not self.t._have_unitcell:
            self.log('No unitcell information')
        else:
            self.log('First frame:')
            self.log('Unitcell angles (deg): %s' % self.t.unitcell_angles[0])
            self.log('Unitcell lengths (nm): %s' % self.t.unitcell_lengths[0])

    def check_topology(self):
        self.section('topology')
        self.log('Number of Atoms:    %d' % ilen(self.t.topology.atoms))
        self.log('Number of Residues: %d' % ilen(self.t.topology.residues))
        self.log('Number of Chains:   %d' % ilen(self.t.topology.chains))
        self.log('Residue names:      %s' % str([r.name for r in self.t.topology.residues]))
        self.log('Unique atom names:  %s' % np.unique([a.name for a in self.t.topology.atoms]))

        # print number of atoms of each element, number of residues of each type?

    def check_bonds(self):
        self.section("Bond Check")
        self.log("Note: PBCs are currently not taken into account during distance check")
        error = False

        for i in range(min(5, self.t.n_frames)):
            #dist = spatial.distance.squareform(spatial.distance.pdist(self.t.xyz[i]))
            for (a, b) in self.t.topology.bonds:
                try:
                    radsum = COVALENT_RADII[a.element.symbol] + COVALENT_RADII[b.element.symbol]
                except KeyError:
                    raise NotImplementedError("I don't have radii information for all of your atoms")

                low, high = self.bond_low * radsum, self.bond_high * radsum
                dist = spatial.distance.euclidean(self.t.xyz[i, a.index], self.t.xyz[i, b.index])
                if not (low < dist < high):
                    error = True
                    self.log('error: atoms %d (%s) and %d (%s) are bonded according '
                            'to the topology but they are a distance of '
                            '%.3f nm apart in frame %d' % (a.index, a.name, b.index, b.name, dist, i))

        self.log("All good.")

    def section(self, title):
        if self._printed_section:
            print()
        print('=== %s: %s ===' % (os.path.basename(self.fn), title))
        self._printed_section = True

    def log(self, msg):
        print(msg)

    def check_positions(self):
        self.section('Positions')
        self.log('Number of frames: %d' % self.t.xyz.shape[0])
        self.log('Number of atoms:  %d' % self.t.xyz.shape[1])
        self.log("Note: PBCs are currently not taken into account during distance check")

        # which atoms have a nonbonded interaction. exluce atoms interacting
        # with themselves or with atoms they're bonded to.
        nonbond_mask = np.logical_not(np.eye(self.t.n_atoms, dtype=np.bool))
        for (a, b) in self.t.topology.bonds:
            nonbond_mask[a.index, b.index] = False
            nonbond_mask[b.index, a.index] = False

        for i in range(min(5, self.t.n_frames)):
            dist = spatial.distance_matrix(self.t.xyz[i], self.t.xyz[i])

            dist = np.ma.array(dist, mask=np.logical_not(nonbond_mask))
            a, b = np.unravel_index(np.ma.argmin(dist), dist.shape)
            names = [q.residue.name + ' ' + q.name for q in self.t.topology.atoms]
            self.log('Frame %d: closest nb dist between '
                     '%d (%s), %d (%s), at d=%.4f nm' % (i, a, names[a], b, names[b], dist[a, b]))

def entry_point():
    args, parser = parse_args()
    main(args, parser)

if __name__ == '__main__':
    entry_point()
