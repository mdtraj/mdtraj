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

import numpy as np
import tempfile
import os
import mdtraj as md
from mdtraj.formats import HDF5TrajectoryFile
from mdtraj.testing import eq
import pytest


try:
    from openmm import unit as units
    HAVE_UNITS = True
except ImportError:
    HAVE_UNITS = False

needs_units = pytest.mark.skipif(not HAVE_UNITS, reason='requires openmm.units')


fd, temp = tempfile.mkstemp(suffix='.h5')
def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.close(fd)
    os.unlink(temp)


def test_write_coordinates():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates)
        assert eq(str(f.root.coordinates.attrs['units']), 'nanometers')


def test_write_coordinates_reshape():
    coordinates = np.random.randn(10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates.reshape(1,10,3))
        assert eq(str(f.root.coordinates.attrs['units']), 'nanometers')


def test_write_multiple():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], np.vstack((coordinates, coordinates)))


def test_write_inconsistent():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)
        # since the first frames we saved didn't contain velocities, we
        # can't save more velocities
        with pytest.raises(ValueError):
            f.write(coordinates, velocities=coordinates)


def test_write_inconsistent_2():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, velocities=coordinates)
        # we're saving a deficient set of data, since before we wrote
        # more information.
        with pytest.raises(ValueError):
            f.write(coordinates)


@needs_units
def test_write_units():
    # openmm.units are automatically converted into MD units for storage on disk
    coordinates = units.Quantity(np.random.randn(4, 10,3), units.angstroms)
    velocities = units.Quantity(np.random.randn(4, 10,3), units.angstroms/units.year)

    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, velocities=velocities)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates.value_in_unit(units.nanometers))
        assert eq(str(f.root.coordinates.attrs['units']), 'nanometers')

        assert eq(f.root.velocities[:], velocities.value_in_unit(units.nanometers/units.picosecond))
        assert eq(str(f.root.velocities.attrs['units']), 'nanometers/picosecond')

def test_write_units2():
    from mdtraj.utils import unit
    coordinates = unit.quantity.Quantity(np.random.randn(4, 10,3),
                    unit.unit_definitions.angstroms)
    velocities = unit.quantity.Quantity(np.random.randn(4, 10,3),
                    unit.unit_definitions.angstroms/unit.unit_definitions.year)

    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, velocities=velocities)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates.value_in_unit(unit.unit_definitions.nanometers))
        assert eq(str(f.root.coordinates.attrs['units']), 'nanometers')

        assert eq(f.root.velocities[:], velocities.value_in_unit(unit.unit_definitions.nanometers/unit.unit_definitions.picosecond))
        assert eq(str(f.root.velocities.attrs['units']), 'nanometers/picosecond')


@needs_units
def test_write_units_mismatch():
    velocities = units.Quantity(np.random.randn(4, 10,3), units.angstroms/units.picosecond)

    with HDF5TrajectoryFile(temp, 'w') as f:
        # if you try to write coordinates that are unitted and not
        # in the correct units, we find that
        with pytest.raises(TypeError):
            f.write(coordinates=velocities)


def test_topology(get_fn):
    top = md.load_pdb(get_fn('native.pdb')).topology

    with HDF5TrajectoryFile(temp, 'w') as f:
        f.topology = top

    with HDF5TrajectoryFile(temp) as f:
        assert f.topology == top


def test_constraints():
    c = np.array([(1,2,3.5)], dtype=np.dtype([('atom1', np.int32), ('atom2', np.int32), ('distance', np.float32)]))

    with HDF5TrajectoryFile(temp, 'w') as f:
        f.constraints = c

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.constraints, c)


def test_constraints2():
    c = np.array([(1,2,3.5)], dtype=np.dtype([('atom1', np.int32), ('atom2', np.int32), ('distance', np.float32)]))

    with HDF5TrajectoryFile(temp, 'w') as f:
        f.constraints = c
        f.constraints = c

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.constraints, c)

def test_read_0():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, alchemicalLambda=np.array([1,2,3,4]))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read()
        assert eq(got.coordinates, coordinates)
        assert eq(got.velocities, None)
        assert eq(got.alchemicalLambda, np.array([1,2,3,4]))


@needs_units
def test_read_1():
    coordinates = units.Quantity(np.random.randn(4, 10,3), units.angstroms)
    velocities = units.Quantity(np.random.randn(4, 10,3), units.angstroms/units.years)


    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, velocities=velocities)

    with HDF5TrajectoryFile(temp) as f:
        got = f.read()
        assert eq(got.coordinates, coordinates.value_in_unit(units.nanometers))
        assert eq(got.velocities, velocities.value_in_unit(units.nanometers/units.picoseconds))


def test_read_slice_0():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, alchemicalLambda=np.array([1,2,3,4]))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(n_frames=2)
        assert eq(got.coordinates, coordinates[:2])
        assert eq(got.velocities, None)
        assert eq(got.alchemicalLambda, np.array([1,2]))


def test_read_slice_1():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(n_frames=2)
        assert eq(got.coordinates, coordinates[:2])
        assert eq(got.velocities, None)

        got = f.read(n_frames=2)
        assert eq(got.coordinates, coordinates[2:])
        assert eq(got.velocities, None)


def test_read_slice_2():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, alchemicalLambda=np.arange(4))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(atom_indices=np.array([0,1]))
        assert eq(got.coordinates, coordinates[:, [0,1], :])
        assert eq(got.alchemicalLambda, np.arange(4))


def test_read_slice_3():
    coordinates = np.random.randn(4, 10,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(coordinates, alchemicalLambda=np.arange(4))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(stride=2, atom_indices=np.array([0,1]))
        assert eq(got.coordinates, coordinates[::2, [0,1], :])
        assert eq(got.alchemicalLambda, np.arange(4)[::2])


def test_do_overwrite():
    with open(temp, 'w') as f:
        f.write('a')

    with HDF5TrajectoryFile(temp, 'w', force_overwrite=True) as f:
        f.write(np.random.randn(10,5,3))


def test_vsite_elements(get_fn):
    #  Test case for issue #265
    pdb_filename = get_fn('GG-tip4pew.pdb')
    trj = md.load(pdb_filename)
    trj.save_hdf5(temp)

    trj2 = md.load(temp, top=pdb_filename)

def test_dont_overwrite():
    with open(temp, 'w') as f:
        f.write('a')
    with pytest.raises(IOError):
        with HDF5TrajectoryFile(temp, 'w', force_overwrite=False) as f:
            f.write(np.random.randn(10,5,3))

def test_attributes():
    constraints = np.zeros(10, dtype=[('atom1', np.int32), ('atom2', np.int32), ('distance', np.float32)])
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.title = 'mytitle'
        f.reference = 'myreference'
        f.forcefield = 'amber99'
        f.randomState = 'sdf'
        f.application = 'openmm'
        f.constraints = constraints

    with HDF5TrajectoryFile(temp) as g:
        eq(g.title, 'mytitle')
        eq(g.reference, 'myreference')
        eq(g.forcefield, 'amber99')
        eq(g.randomState, 'sdf')
        eq(g.application, 'openmm')
        eq(g.constraints, constraints)

def test_append():
    x1 = np.random.randn(10,5,3)
    x2 = np.random.randn(8,5,3)
    with HDF5TrajectoryFile(temp, 'w') as f:
        f.write(x1)
    with HDF5TrajectoryFile(temp, 'a') as f:
        f.write(x2)

    with HDF5TrajectoryFile(temp) as f:
        eq(f.root.coordinates[:], np.concatenate((x1,x2)))
