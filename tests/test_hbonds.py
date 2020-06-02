##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2017 Stanford University and the Authors
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


import os
import subprocess
from distutils.spawn import find_executable

import mdtraj as md
import numpy as np
import pytest
import scipy.sparse
from mdtraj.testing import eq

DSSP_MSG = "This test requires mkdssp to be installed, from http://swift.cmbi.ru.nl/gv/dssp/"
needs_dssp = pytest.mark.skipif(not find_executable('mkdssp'), reason=DSSP_MSG)


def test_hbonds(get_fn):
    t = md.load(get_fn('2EQQ.pdb'))
    ours = md.geometry.hbond.kabsch_sander(t)


@needs_dssp
def test_hbonds_against_dssp(get_fn, tmpdir):
    t = md.load(get_fn('2EQQ.pdb'))[0]
    pdb = os.path.join(tmpdir, 'f.pdb')
    dssp = os.path.join(tmpdir, 'f.pdb.dssp')
    t.save(pdb)

    cmd = ['mkdssp', '-i', pdb, '-o', dssp]
    subprocess.check_output(' '.join(cmd), shell=True)
    energy = scipy.sparse.lil_matrix((t.n_residues, t.n_residues))

    # read the dssp N-H-->O column from the output file
    with open(dssp) as f:
        # skip the lines until the description of each residue's hbonds
        while not f.readline().startswith('  #  RESIDUE AA STRUCTURE'):
            continue

        for i, line in enumerate(f):
            line = line.rstrip()
            offset0, e0 = map(float, line[39:50].split(','))
            offset1, e1 = map(float, line[61:72].split(','))
            if e0 <= -0.5:
                energy[int(i + offset0), i] = e0
            if e1 <= -0.5:
                energy[int(i + offset1), i] = e1

    dssp = energy.todense()
    ours = md.geometry.hbond.kabsch_sander(t)[0].todense()

    # There is tricky issues with the rounding right at the -0.5 cutoff,
    # so lets just check for equality with DSSP at -0.6 or less
    eq((dssp < -0.6), (ours < -0.6))
    eq(dssp[dssp < -0.6], ours[ours < -0.6], decimal=1)


def test_baker_hubbard_0(get_fn):
    t = md.load(get_fn('2EQQ.pdb'))

    # print('to view the hbonds defined in 2EQQ by baker_hubbard()')
    # print('put these commands into pymol on top of the pdb:\n')
    # for e in md.geometry.hbond.baker_hubbard(t):
    #     print('distance RANK %d, RANK %d' % (e[1], e[2]))

    # these are the results produced by the algorithm on this protein as
    # of 11/26/13. This unit test basically just ensures that the method
    # runs and produces the same results it did then. It's no guarentee that
    # these are the "TRUE" hydrogen bonds in this system.
    ref = np.array([[0, 10, 8], [0, 11, 7], [69, 73, 54], [76, 82, 65],
                    [119, 131, 89], [140, 148, 265], [166, 177, 122],
                    [181, 188, 231], [209, 217, 215], [221, 225, 184],
                    [228, 239, 186], [235, 247, 216], [262, 271, 143],
                    [298, 305, 115], [186, 191, 215], [413, 419, 392]])
    eq(ref, md.geometry.hbond.baker_hubbard(t))


def test_baker_hubbard_1(get_fn):
    # no hydrogens in this file -> no hydrogen bonds
    t = md.load(get_fn('1bpi.pdb'))
    eq(np.zeros((0, 3), dtype=int), md.baker_hubbard(t))


def test_baker_hubbard_2(get_fn):
    t = md.load(get_fn('1vii_sustiva_water.pdb'))
    triplets = md.baker_hubbard(t)
    N = 1000
    rows = triplets[:, 0] * N * N + triplets[:, 1] * N + triplets[:, 2]
    # ensure that there aren't any repeat rows
    eq(len(np.unique(rows)), len(rows))


def test_baker_hubbard_3(get_fn):
    #different distance cutoffs->different hydrogen bonds found
    t = md.load(get_fn('2waters_baker_hubbard.pdb'))
    eq(np.zeros((0,3), dtype=int), md.baker_hubbard(t, exclude_water=False, distance_cutoff = 0.18))
    eq(np.zeros((0,3), dtype=int), md.baker_hubbard(t, exclude_water=False, angle_cutoff=180))
    eq(np.array([[0, 1, 3]]), md.baker_hubbard(t, exclude_water=False))


def test_wernet_nilsson_0(get_fn):
    # no hydrogens in this file -> no hydrogen bonds
    t0 = md.load(get_fn('1bpi.pdb'))
    assert len(md.wernet_nilsson(t0)) == len(t0)
    eq(np.zeros((0, 3), dtype=int), md.wernet_nilsson(t0)[0])


def test_wernet_nilsson_1(get_fn):
    # one of these files has PBCs and the other doesnt
    for fn in ['2EQQ.pdb', '4ZUO.pdb']:
        t = md.load(get_fn(fn))

        result = md.wernet_nilsson(t)
        assert len(result) == len(t)
        assert isinstance(result, list)
        assert all(isinstance(elem, np.ndarray) for elem in result)
        assert all(elem.shape[1] == 3 for elem in result)
        for frame, hbonds in enumerate(result):
            for d_i, h_i, a_i in hbonds:
                assert t.topology.atom(d_i).element.symbol in ['O', 'N']
                assert t.topology.atom(h_i).element.symbol == 'H'
                assert t.topology.atom(a_i).element.symbol in ['O', 'N']
            # assert that the donor-acceptor distance is less than 0.5 nm, just
            # to make sure the criterion is giving back totally implausible stuff
            if len(hbonds) > 0:
                assert np.all(md.compute_distances(t[frame], hbonds[:, [0, 2]]) < 0.5)
