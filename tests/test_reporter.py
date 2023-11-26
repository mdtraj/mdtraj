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
import shutil

import numpy as np
import pytest

import mdtraj as md
from mdtraj.formats import HDF5TrajectoryFile, NetCDFTrajectoryFile
from mdtraj.reporters import HDF5Reporter, NetCDFReporter, DCDReporter, XTCReporter
from mdtraj.testing import eq

try:
    from openmm.unit import nanometers, kelvin, picoseconds, femtoseconds
    from openmm import LangevinIntegrator, Platform
    from openmm.app import PDBFile, ForceField, Simulation, CutoffNonPeriodic, CutoffPeriodic, HBonds, CheckpointReporter

    HAVE_OPENMM = True
except ImportError:
    HAVE_OPENMM = False

# special pytest global to mark all tests in this module
pytestmark = pytest.mark.skipif(not HAVE_OPENMM, reason='test_reporter.py needs OpenMM.')


def test_reporter(tmpdir, get_fn):
    pdb = PDBFile(get_fn('native.pdb'))
    forcefield = ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    # NO PERIODIC BOUNDARY CONDITIONS
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                                     nonbondedCutoff=1.0 * nanometers, constraints=HBonds, rigidWater=True)
    integrator = LangevinIntegrator(300 * kelvin, 1.0 / picoseconds, 2.0 * femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = Platform.getPlatformByName('Reference')
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    simulation.context.setVelocitiesToTemperature(300 * kelvin)

    tmpdir = str(tmpdir)
    hdf5file = os.path.join(tmpdir, 'traj.h5')
    ncfile = os.path.join(tmpdir, 'traj.nc')
    dcdfile = os.path.join(tmpdir, 'traj.dcd')
    xtcfile = os.path.join(tmpdir, 'traj.xtc')

    reporter = HDF5Reporter(hdf5file, 2, coordinates=True, time=True,
                            cell=True, potentialEnergy=True, kineticEnergy=True, temperature=True,
                            velocities=True)
    reporter2 = NetCDFReporter(ncfile, 2, coordinates=True, time=True, cell=True)
    reporter3 = DCDReporter(dcdfile, 2)
    reporter4 = XTCReporter(xtcfile, 2)

    simulation.reporters.append(reporter)
    simulation.reporters.append(reporter2)
    simulation.reporters.append(reporter3)
    simulation.reporters.append(reporter4)
    simulation.step(100)

    reporter.close()
    reporter2.close()
    reporter3.close()
    reporter4.close()

    with HDF5TrajectoryFile(hdf5file) as f:
        got = f.read()
        eq(got.temperature.shape, (50,))
        eq(got.potentialEnergy.shape, (50,))
        eq(got.kineticEnergy.shape, (50,))
        eq(got.coordinates.shape, (50, 22, 3))
        eq(got.velocities.shape, (50, 22, 3))
        eq(got.cell_lengths.shape, (50, 3))
        eq(got.cell_angles.shape, (50, 3))
        eq(got.time, 0.002 * 2 * (1 + np.arange(50)))
        assert f.topology == md.load(get_fn('native.pdb')).top

    with NetCDFTrajectoryFile(ncfile) as f:
        xyz, time, cell_lengths, cell_angles = f.read()
        eq(cell_lengths.shape, (50, 3))
        eq(cell_angles.shape, (50, 3))
        eq(time, 0.002 * 2 * (1 + np.arange(50)))

    hdf5_traj = md.load(hdf5file)
    dcd_traj = md.load(dcdfile, top=get_fn('native.pdb'))
    netcdf_traj = md.load(ncfile, top=get_fn('native.pdb'))
    xtc_traj = md.load(xtcfile, top=get_fn('native.pdb'))

    # we don't have to convert units here, because md.load already
    # handles that
    eq(hdf5_traj.xyz, netcdf_traj.xyz)
    eq(hdf5_traj.unitcell_vectors, netcdf_traj.unitcell_vectors)
    eq(hdf5_traj.time, netcdf_traj.time)
    eq(xtc_traj.time, netcdf_traj.time)

    eq(dcd_traj.xyz, hdf5_traj.xyz)
    eq(xtc_traj.xyz, dcd_traj.xyz, decimal=3)
    eq(dcd_traj.unitcell_vectors, hdf5_traj.unitcell_vectors)


def test_reporter_subset(tmpdir, get_fn):
    pdb = PDBFile(get_fn('native2.pdb'))
    pdb.topology.setUnitCellDimensions([2, 2, 2])
    forcefield = ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffPeriodic,
                                     nonbondedCutoff=1 * nanometers, constraints=HBonds, rigidWater=True)
    integrator = LangevinIntegrator(300 * kelvin, 1.0 / picoseconds, 2.0 * femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = Platform.getPlatformByName('Reference')
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    simulation.context.setVelocitiesToTemperature(300 * kelvin)

    tmpdir = str(tmpdir)
    hdf5file = os.path.join(tmpdir, 'traj.h5')
    ncfile = os.path.join(tmpdir, 'traj.nc')
    dcdfile = os.path.join(tmpdir, 'traj.dcd')
    xtcfile = os.path.join(tmpdir, 'traj.xtc')

    atomSubset = [0, 1, 2, 4, 5]

    reporter = HDF5Reporter(hdf5file, 2, coordinates=True, time=True,
                            cell=True, potentialEnergy=True, kineticEnergy=True, temperature=True,
                            velocities=True, atomSubset=atomSubset)
    reporter2 = NetCDFReporter(ncfile, 2, coordinates=True, time=True,
                               cell=True, atomSubset=atomSubset)
    reporter3 = DCDReporter(dcdfile, 2, atomSubset=atomSubset)
    reporter4 = XTCReporter(xtcfile, 2, atomSubset=atomSubset)

    simulation.reporters.append(reporter)
    simulation.reporters.append(reporter2)
    simulation.reporters.append(reporter3)
    simulation.reporters.append(reporter4)
    simulation.step(100)

    reporter.close()
    reporter2.close()
    reporter3.close()
    reporter4.close()

    t = md.load(get_fn('native.pdb'))
    t.restrict_atoms(atomSubset)

    with HDF5TrajectoryFile(hdf5file) as f:
        got = f.read()
        eq(got.temperature.shape, (50,))
        eq(got.potentialEnergy.shape, (50,))
        eq(got.kineticEnergy.shape, (50,))
        eq(got.coordinates.shape, (50, len(atomSubset), 3))
        eq(got.velocities.shape, (50, len(atomSubset), 3))
        eq(got.cell_lengths, 2 * np.ones((50, 3)))
        eq(got.cell_angles, 90 * np.ones((50, 3)))
        eq(got.time, 0.002 * 2 * (1 + np.arange(50)))
        assert f.topology == md.load(get_fn('native.pdb'), atom_indices=atomSubset).topology

    with NetCDFTrajectoryFile(ncfile) as f:
        xyz, time, cell_lengths, cell_angles = f.read()
        eq(cell_lengths, 20 * np.ones((50, 3)))
        eq(cell_angles, 90 * np.ones((50, 3)))
        eq(time, 0.002 * 2 * (1 + np.arange(50)))
        eq(xyz.shape, (50, len(atomSubset), 3))

    hdf5_traj = md.load(hdf5file)
    dcd_traj = md.load(dcdfile, top=hdf5_traj)
    netcdf_traj = md.load(ncfile, top=hdf5_traj)
    xtc_traj = md.load(xtcfile, top=hdf5_traj)

    # we don't have to convert units here, because md.load already handles that
    eq(hdf5_traj.xyz, netcdf_traj.xyz)
    eq(hdf5_traj.unitcell_vectors, netcdf_traj.unitcell_vectors)
    eq(hdf5_traj.time, netcdf_traj.time)
    eq(xtc_traj.time, netcdf_traj.time)

    eq(dcd_traj.xyz, hdf5_traj.xyz)
    eq(xtc_traj.xyz, hdf5_traj.xyz)
    eq(dcd_traj.unitcell_vectors, hdf5_traj.unitcell_vectors)


def test_xtc_reporter_append(tmpdir, get_fn):
    pdb = PDBFile(get_fn('native.pdb'))
    forcefield = ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    # NO PERIODIC BOUNDARY CONDITIONS
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                                     nonbondedCutoff=1.0 * nanometers, constraints=HBonds, rigidWater=True)
    integrator = LangevinIntegrator(300 * kelvin, 1.0 / picoseconds, 2.0 * femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = Platform.getPlatformByName('Reference')
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    simulation.context.setVelocitiesToTemperature(300 * kelvin)

    tmpdir = str(tmpdir)
    xtcfile = os.path.join(tmpdir, 'traj.xtc')
    xtcfile_cp = os.path.join(tmpdir, 'traj_cp.xtc')
    checkpoint = os.path.join(tmpdir, 'checkpoint.chk')
    reporter = XTCReporter(xtcfile, 2)
    simulation.reporters.append(reporter)
    simulation.reporters.append(CheckpointReporter(checkpoint, 10))
    simulation.step(10)
    reporter.close()
    shutil.copyfile(xtcfile, xtcfile_cp)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                                     nonbondedCutoff=1.0 * nanometers, constraints=HBonds, rigidWater=True)
    integrator = LangevinIntegrator(300 * kelvin, 1.0 / picoseconds, 2.0 * femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = Platform.getPlatformByName('Reference')
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.loadCheckpoint(checkpoint)
    reporter = XTCReporter(xtcfile, 2, append=True)
    simulation.reporters.append(reporter)
    simulation.step(10)
    reporter.close()
    xtc_traj = md.load(xtcfile, top=get_fn('native.pdb'))
    xtc_traj_cp = md.load(xtcfile_cp, top=get_fn('native.pdb'))
    eq(xtc_traj.xyz[:5], xtc_traj_cp.xyz)
    eq(xtc_traj.n_frames, 10)
    eq(xtc_traj_cp.n_frames, 5)
    eq(xtc_traj.time[:5], xtc_traj_cp.time)
