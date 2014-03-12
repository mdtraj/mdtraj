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
from mdtraj.core.trajectory import _parse_topology

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
        in any supported format.''')
    # parser.add_argument('-n', '--noload', action='store_true', help='''Do not load the coordinate data from the trajectory, for example if the trajectory is too large to load into memory. Only a limited number of checks will be done.''')
    parser.add_argument('-t', '--topology', type=str, help='''Topology for the system (.prmtop/.pdb)''')
    parser.add_argument('--bond-low', type=float, help='''Minimum fraction of sum of covalent radii for bonded atoms. Default=0.4''', default=0.4)
    parser.add_argument('--bond-high', type=float, help='''Maximum fraction of sum of covalent radii for bonded atoms. Default=1.2''', default=1.2)
    parser.add_argument('--rmsd-tolerance', type=float, help='''Maximum tolerance for percent change in RMSD. Default=100.0''', default=100.0)
    return parser.parse_args(), parser


def main(args, parser):
    inspector = Inspector(args.bond_low, args.bond_high, args.rmsd_tolerance)

    topology_files = [fn for fn in args.files if os.path.splitext(fn)[1] in ['.pdb', '.prmtop', '.h5']]
    if args.topology is None and len(topology_files) > 0:
        args.topology = topology_files[0]

    inspector.load_topology(args.topology)

    for f in args.files:
        if not os.path.exists(f):
            parser.error("File '%s' does not exist" % f)
        if not os.path.isfile(f):
            parser.error("File '%s' is not a file" % f)

        inspector.load_trajectory(f)


class Inspector(object):
    def __init__(self, bond_low, bond_high, rmsd_tolerance):
        self._printed_section = False
        self.t = None
        self.fn = None

        self.bond_low = bond_low
        self.bond_high = bond_high
        self.rmsd_tolerance = rmsd_tolerance

    def load_topology(self, fn):
        self.fn = fn
        self.topology = _parse_topology(fn)
        self.check_topology()

    def load_trajectory(self, fn):
        self.fn = fn

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.t = md.load(fn, top=self.topology)

        self.check_shape()
        self.check_unitcell()
        self.check_bonds()
        self.check_nonbonded()
        self.check_imaging()

    def check_shape(self):
        self.section('Shape')
        self.log('Number of frames: %d' % self.t.xyz.shape[0])
        self.log('Number of atoms:  %d' % self.t.xyz.shape[1])

    def check_unitcell(self):
        self.section('Unitcell')
        if not self.t._have_unitcell:
            self.log('No unitcell information')
        else:
            self.log('First frame:')
            self.log('Unitcell angles (deg): %s' % self.t.unitcell_angles[0])
            self.log('Unitcell lengths (nm): %s' % self.t.unitcell_lengths[0])

    def check_topology(self):
        self.section('topology')
        self.log('Number of Atoms:    %d' % ilen(self.topology.atoms))
        self.log('Number of Residues: %d' % ilen(self.topology.residues))
        self.log('Number of Chains:   %d' % ilen(self.topology.chains))
        self.log('Residues:           %s' % ', '.join(['%s (%d atoms)' % (r, ilen(r.atoms)) for r in self.topology.residues]))
        self.log('Unique atom names:  %s' % ', '.join(np.unique([a.name for a in self.topology.atoms])))

    def check_bonds(self):
        self.section("Bond Check (without PBCs)")

        radii = []
        pairs = []
        for (a, b) in self.topology.bonds:
            try:
                radsum = COVALENT_RADII[a.element.symbol] + COVALENT_RADII[b.element.symbol]
            except KeyError:
                raise NotImplementedError("I don't have radii information for all of your atoms")
            radii.append(radsum)
            pairs.append((a.index, b.index))

        radii = np.array(radii)
        pairs = np.array(pairs)

        distances = md.compute_distances(self.t, pairs, periodic=False)
        low, high = self.bond_low * radii, self.bond_high * radii
        extreme = np.logical_or(distances < low, distances > high)

        if np.any(extreme):
            frames, bonds = np.nonzero(extreme)
            frame, bond = frames[0], bonds[0]
            a1 = self.topology.atom(pairs[bond][0])
            a2 = self.topology.atom(pairs[bond][0])

            self.log('error: atoms (%s) and (%s) are bonded according to the topology ' % (a1, a2))
            self.log('but they are a distance of %.3f nm apart in frame %d' % (distances[frame, bond], frame))
        else:
            self.log("All good.")

    def section(self, title):
        if self._printed_section:
            print()
        print('=== %s: %s ===' % (os.path.basename(self.fn), title))
        self._printed_section = True

    def log(self, msg):
        print(msg)

    def check_nonbonded(self):
        self.section('Nonbonded Check (without PBCs)')

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

    def check_imaging(self):
        if self.t.n_frames > 2:
            self.section('Imaging')

            r = md.rmsd(target=self.t, reference=self.t[0])
            percent_change = np.divide(np.abs(r[2:]-r[1:-1]), r[1:-1]) * 100.0
            if np.any(percent_change > self.rmsd_tolerance):
                self.log('Potential imaging issue: %s' % self.fn)
            else:
                self.log('No imaging issue detected')


def entry_point():
    args, parser = parse_args()
    main(args, parser)


if __name__ == '__main__':
    entry_point()
