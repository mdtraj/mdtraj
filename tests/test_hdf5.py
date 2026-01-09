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

import os
import tempfile

import numpy as np
import pytest

import mdtraj as md
from mdtraj.core.topology import Topology
from mdtraj.formats import HDF5TrajectoryFile
from mdtraj.testing import eq

try:
    from openmm import unit as units

    HAVE_UNITS = True
except ImportError:
    HAVE_UNITS = False

needs_units = pytest.mark.skipif(not HAVE_UNITS, reason="requires openmm.units")


fd, temp = tempfile.mkstemp(suffix=".h5")
rng = np.random.default_rng()


def teardown_module(module):
    """remove the temporary file created by tests in this file
    this gets automatically called by pytest"""
    os.close(fd)
    os.unlink(temp)


def test_write_coordinates():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates)
        assert eq(str(f.root.coordinates.attrs["units"]), "nanometers")


def test_write_coordinates_reshape():
    coordinates = rng.standard_normal((10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates.reshape(1, 10, 3))
        assert eq(str(f.root.coordinates.attrs["units"]), "nanometers")


def test_write_multiple():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates)
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], np.vstack((coordinates, coordinates)))


def test_write_inconsistent():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates)
        # since the first frames we saved didn't contain velocities, we
        # can't save more velocities
        with pytest.raises(ValueError):
            f.write(coordinates, velocities=coordinates)


def test_write_inconsistent_2():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, velocities=coordinates)
        # we're saving a deficient set of data, since before we wrote
        # more information.
        with pytest.raises(ValueError):
            f.write(coordinates)


@needs_units
def test_write_units():
    # openmm.units are automatically converted into MD units for storage on disk
    coordinates = units.Quantity(rng.standard_normal((4, 10, 3)), units.angstroms)
    velocities = units.Quantity(rng.standard_normal((4, 10, 3)), units.angstroms / units.year)

    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, velocities=velocities)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.root.coordinates[:], coordinates.value_in_unit(units.nanometers))
        assert eq(str(f.root.coordinates.attrs["units"]), "nanometers")

        assert eq(
            f.root.velocities[:],
            velocities.value_in_unit(units.nanometers / units.picosecond),
        )
        assert eq(str(f.root.velocities.attrs["units"]), "nanometers/picosecond")


def test_write_units2():
    from mdtraj.utils import unit

    coordinates = unit.quantity.Quantity(
        rng.standard_normal((4, 10, 3)),
        unit.unit_definitions.angstroms,
    )
    velocities = unit.quantity.Quantity(
        rng.standard_normal((4, 10, 3)),
        unit.unit_definitions.angstroms / unit.unit_definitions.year,
    )

    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, velocities=velocities)

    with HDF5TrajectoryFile(temp) as f:
        assert eq(
            f.root.coordinates[:],
            coordinates.value_in_unit(unit.unit_definitions.nanometers),
        )
        assert eq(str(f.root.coordinates.attrs["units"]), "nanometers")

        assert eq(
            f.root.velocities[:],
            velocities.value_in_unit(
                unit.unit_definitions.nanometers / unit.unit_definitions.picosecond,
            ),
        )
        assert eq(str(f.root.velocities.attrs["units"]), "nanometers/picosecond")


@needs_units
def test_write_units_mismatch():
    velocities = units.Quantity(
        rng.standard_normal((4, 10, 3)),
        units.angstroms / units.picosecond,
    )

    with HDF5TrajectoryFile(temp, "w") as f:
        # if you try to write coordinates that are unitted and not
        # in the correct units, we find that
        with pytest.raises(TypeError):
            f.write(coordinates=velocities)


def test_topology(get_fn):
    top = md.load_pdb(get_fn("native.pdb")).topology

    with HDF5TrajectoryFile(temp, "w") as f:
        f.topology = top

    with HDF5TrajectoryFile(temp) as f:
        assert f.topology == top


def test_constraints():
    c = np.array(
        [(1, 2, 3.5)],
        dtype=np.dtype(
            [("atom1", np.int32), ("atom2", np.int32), ("distance", np.float32)],
        ),
    )

    with HDF5TrajectoryFile(temp, "w") as f:
        f.constraints = c

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.constraints, c)


def test_constraints2():
    c = np.array(
        [(1, 2, 3.5)],
        dtype=np.dtype(
            [("atom1", np.int32), ("atom2", np.int32), ("distance", np.float32)],
        ),
    )

    with HDF5TrajectoryFile(temp, "w") as f:
        f.constraints = c
        f.constraints = c

    with HDF5TrajectoryFile(temp) as f:
        assert eq(f.constraints, c)


def test_read_0():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, alchemicalLambda=np.array([1, 2, 3, 4]))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read()
        assert eq(got.coordinates, coordinates)
        assert eq(got.velocities, None)
        assert eq(got.alchemicalLambda, np.array([1, 2, 3, 4]))


@needs_units
def test_read_1():
    coordinates = units.Quantity(rng.standard_normal((4, 10, 3)), units.angstroms)
    velocities = units.Quantity(
        rng.standard_normal((4, 10, 3)),
        units.angstroms / units.years,
    )

    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, velocities=velocities)

    with HDF5TrajectoryFile(temp) as f:
        got = f.read()
        assert eq(got.coordinates, coordinates.value_in_unit(units.nanometers))
        assert eq(
            got.velocities,
            velocities.value_in_unit(units.nanometers / units.picoseconds),
        )


def test_read_slice_0():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, alchemicalLambda=np.array([1, 2, 3, 4]))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(n_frames=2)
        assert eq(got.coordinates, coordinates[:2])
        assert eq(got.velocities, None)
        assert eq(got.alchemicalLambda, np.array([1, 2]))


def test_read_slice_1():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates)

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(n_frames=2)
        assert eq(got.coordinates, coordinates[:2])
        assert eq(got.velocities, None)

        got = f.read(n_frames=2)
        assert eq(got.coordinates, coordinates[2:])
        assert eq(got.velocities, None)


def test_read_slice_2():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, alchemicalLambda=np.arange(4))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(atom_indices=np.array([0, 1]))
        assert eq(got.coordinates, coordinates[:, [0, 1], :])
        assert eq(got.alchemicalLambda, np.arange(4))


def test_read_slice_3():
    coordinates = rng.standard_normal((4, 10, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(coordinates, alchemicalLambda=np.arange(4))

    with HDF5TrajectoryFile(temp) as f:
        got = f.read(stride=2, atom_indices=np.array([0, 1]))
        assert eq(got.coordinates, coordinates[::2, [0, 1], :])
        assert eq(got.alchemicalLambda, np.arange(4)[::2])


def test_do_overwrite():
    with open(temp, "w") as f:
        f.write("a")

    with HDF5TrajectoryFile(temp, "w", force_overwrite=True) as f:
        f.write(rng.standard_normal((10, 5, 3)))


def test_vsite_elements(get_fn):
    #  Test case for issue #265
    pdb_filename = get_fn("GG-tip4pew.pdb")
    trj = md.load(pdb_filename)
    trj.save_hdf5(temp)

    md.load(temp, top=pdb_filename)


def test_dont_overwrite():
    with open(temp, "w") as f:
        f.write("a")
    with pytest.raises(IOError):
        with HDF5TrajectoryFile(temp, "w", force_overwrite=False) as f:
            f.write(rng.standard_normal((10, 5, 3)))


def test_attributes():
    constraints = np.zeros(
        10,
        dtype=[("atom1", np.int32), ("atom2", np.int32), ("distance", np.float32)],
    )
    with HDF5TrajectoryFile(temp, "w") as f:
        f.title = "mytitle"
        f.reference = "myreference"
        f.forcefield = "amber99"
        f.randomState = "sdf"
        f.application = "openmm"
        f.constraints = constraints

    with HDF5TrajectoryFile(temp) as g:
        eq(g.title, "mytitle")
        eq(g.reference, "myreference")
        eq(g.forcefield, "amber99")
        eq(g.randomState, "sdf")
        eq(g.application, "openmm")
        eq(g.constraints, constraints)


def test_append():
    x1 = rng.standard_normal((10, 5, 3))
    x2 = rng.standard_normal((8, 5, 3))
    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(x1)
    with HDF5TrajectoryFile(temp, "a") as f:
        f.write(x2)

    with HDF5TrajectoryFile(temp) as f:
        eq(f.root.coordinates[:], np.concatenate((x1, x2)))


def test_topology_None(h5traj):
    with HDF5TrajectoryFile(temp, "w") as f:
        f.topology = None

        assert f.topology is None, "The None argument did not pass through the topology setter properly"

    _, in_fn, _ = h5traj
    with HDF5TrajectoryFile(in_fn, "a") as f:
        # The file previously has a topology
        assert isinstance(f.topology, Topology), "The test HDF5 File does not contain a topology"
        f.topology = None

        # The topology should now be overwritten as None now
        assert f.topology is None, "The topology of the HDF5 file was not deleted"


def test_hdf5_bond_metadata(get_fn):
    # test that bond metadata is preserved when writing
    # and reading HDF5. Added in PR 2101

    traj = md.load(get_fn("imatinib.pdb"), bond_orders=True)

    with HDF5TrajectoryFile(temp, "w", bond_metadata=True) as f:
        f.write(traj.xyz)
        f.topology = traj.topology

    # Check that bond metadata was saved correctly
    with HDF5TrajectoryFile(temp) as f:
        top = f.topology
        for bond1, bond2 in zip(traj.topology.bonds, top.bonds):
            assert bond1 == bond2


def test_hdf5_formal_charge(get_fn):
    # test that atom formal charge is preserved when writing
    # and reading HDF5. Added in PR 2101

    traj = md.load(get_fn("1ply_charge.pdb"))

    with HDF5TrajectoryFile(temp, "w") as f:
        f.write(traj.xyz)
        f.topology = traj.topology

    # Check that formal charges were saved correctly
    with HDF5TrajectoryFile(temp) as f:
        top = f.topology
        for atom1, atom2 in zip(traj.topology.atoms, top.atoms):
            assert atom1.formal_charge == atom2.formal_charge


def test_hdf5_roundtrip(h5traj_full_metadata):
    # test that formal_charge and bond order/type is preserved
    # on topology when saving to and reading from HDF5.

    traj, _, tmp_dir = h5traj_full_metadata
    out_fn = os.path.join(tmp_dir, "roundtrip.h5")

    traj.save(out_fn)
    traj2 = md.load(out_fn)

    assert eq(traj.xyz, traj2.xyz)
    assert traj.topology == traj2.topology
    assert eq(traj.time, traj2.time)
    assert eq(traj.unitcell_lengths, traj2.unitcell_lengths)
    assert eq(traj.unitcell_angles, traj2.unitcell_angles)
    for atom1, atom2 in zip(traj.topology.atoms, traj2.topology.atoms):
        assert atom1.formal_charge == atom2.formal_charge
    for bond1, bond2 in zip(traj.topology.bonds, traj2.topology.bonds):
        assert bond1.order == bond2.order
        assert bond1.type == bond2.type
