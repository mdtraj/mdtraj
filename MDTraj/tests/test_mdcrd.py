# Copyright 2012 mdtraj developers
#
# This file is part of mdtraj
#
# mdtraj is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mdtraj is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mdtraj. If not, see http://www.gnu.org/licenses/.

import tempfile, os
import numpy as np
import mdtraj as md
from mdtraj import MDCRDTrajectoryFile, mdcrd
from mdtraj.testing import get_fn, eq, DocStringFormatTester, raises
TestDocstrings = DocStringFormatTester(mdcrd, error_on_none=True)

temp = tempfile.mkstemp(suffix='.mdcrd')[1]
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by nose"""
    os.unlink(temp)


def test_read_0():
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrd'), n_atoms=22) as f:
        xyz, box = f.read()
        assert box is None
    with MDCRDTrajectoryFile(get_fn('frame0.mdcrdbox'), n_atoms=22) as f:
        xyz2, box = f.read()
        eq(box.shape, (1002, 3))
        eq(xyz, xyz2)


def test_read_write_0():
    xyz = 10*np.random.randn(100, 11, 3)
    with MDCRDTrajectoryFile(temp, mode='w') as f:
        f.write(xyz)
    with MDCRDTrajectoryFile(temp, n_atoms=11) as f:
        xyz2, _ = f.read()

    eq(_, None)
    eq(xyz, xyz2, decimal=3)


def test_read_write_1():
    xyz = 10*np.random.randn(100, 11, 3)
    box = np.random.randn(100,3)
    with MDCRDTrajectoryFile(temp, mode='w') as f:
        f.write(xyz, box)
    with MDCRDTrajectoryFile(temp, n_atoms=11) as f:
        xyz2, box2 = f.read()

    eq(box, box2, decimal=3)
    eq(xyz, xyz2, decimal=3)


def test_read_write_2():
    pdb = md.load(get_fn('1bpi.pdb'))
    pdb.save(temp)

    t = md.load(temp, top=pdb.topology)

    eq(t.xyz, pdb.xyz)
    eq(t.unitcell_vectors, pdb.unitcell_vectors)
