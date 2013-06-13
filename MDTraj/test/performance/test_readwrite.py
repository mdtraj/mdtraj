import os
import tempfile
import numpy as np
from noseperf.testcases import PerformanceTest

from mdtraj import xtc, trr, netcdf, hdf5
from mdtraj import DCDTrajectoryFile, BINPOSTrajectoryFile

######################################################
# Base class: handles setting up a temp file and array
######################################################


class WithTemp(PerformanceTest):
    n_frames = 10000
    n_atoms = 100

    def setUp(self):
        self.xyz = np.random.randn(self.n_frames, self.n_atoms, 3).astype(np.float32)
        self.fn = tempfile.mkstemp()[1]

    def tearDown(self):
        os.unlink(self.fn)

########################################
# Tests
########################################


class TestXTCWrite(WithTemp):
    def test(self):
        "Test the write speed of the XTC code (10000 frames, 100 atoms)"
        xtc.write(self.fn, xyz=self.xyz, force_overwrite=True)


class TestXTCRead(WithTemp):
    def setUp(self):
        super(TestXTCRead, self).setUp()
        xtc.write(self.fn, xyz=self.xyz, force_overwrite=True)

    def test(self):
        "Test the read speed of the XTC code (10000 frames, 100 atoms)"
        xtc.read(self.fn)


class TestDCDWrite(WithTemp):
    def test(self):
        "Test the write speed of the DCD code (10000 frames, 100 atoms)"
        with DCDTrajectoryFile(self.fn, 'w') as f:
            f.write(xyz=self.xyz)


class TestDCDRead(WithTemp):
    def setUp(self):
        super(TestDCDRead, self).setUp()
        with DCDTrajectoryFile(self.fn, 'w') as f:
            f.write(xyz=self.xyz)

    def test(self):
        "Test the read speed of the DCD code (10000 frames, 100 atoms)"
        with DCDTrajectoryFile(self.fn) as f:
            f.read()


class TestBINPOSWrite(WithTemp):
    def test(self):
        "Test the write speed of the BINPOS code (10000 frames, 100 atoms)"
        with BINPOSTrajectoryFile(self.fn, 'w', force_overwrite=True) as f:
            f.write(self.xyz)


class TestBINPOSRead(WithTemp):
    def setUp(self):
        super(TestBINPOSRead, self).setUp()
        with BINPOSTrajectoryFile(self.fn, 'w', force_overwrite=True) as f:
            f.write(self.xyz)

    def test(self):
        "Test the read speed of the BINPOS code (10000 frames, 100 atoms)"
        with BINPOSTrajectoryFile(self.fn) as f:
            xyz = f.read()


class TestTRRWrite(WithTemp):
    def test(self):
        "Test the write speed of the TRR code (10000 frames, 100 atoms)"
        trr.write(self.fn, xyz=self.xyz, force_overwrite=True)


class TestTRRRead(WithTemp):
    def setUp(self):
        super(TestTRRRead, self).setUp()
        trr.write(self.fn, xyz=self.xyz, force_overwrite=True)

    def test(self):
        "Test the read speed of the TRR code (10000 frames, 100 atoms)"
        trr.read(self.fn)


class TestNetCDFWrite(WithTemp):
    def test(self):
        "Test the write speed of the NetCDF code (10000 frames, 100 atoms)"
        with netcdf.NetCDFTrajectoryFile(self.fn, 'w', force_overwrite=True) as f:
            f.write(self.xyz)


class TestNetCDFRead(WithTemp):
    def setUp(self):
        super(TestNetCDFRead, self).setUp()
        with netcdf.NetCDFTrajectoryFile(self.fn, 'w', force_overwrite=True) as f:
            f.write(self.xyz)

    def test(self):
        "Test the read speed of the NetCDF code (10000 frames, 100 atoms)"
        with netcdf.NetCDFTrajectoryFile(self.fn) as f:
            f.read()


class TestHDF5Write(WithTemp):
    def test(self):
        "Test the write speed of the hdf5 code (10000 frames, 100 atoms)"
        with hdf5.HDF5TrajectoryFile(self.fn, 'w', force_overwrite=True) as f:
            f.write(self.xyz)


class TestHDF5Read(WithTemp):
    def setUp(self):
        super(TestHDF5Read, self).setUp()
        with hdf5.HDF5TrajectoryFile(self.fn, 'w', force_overwrite=True) as f:
            f.write(self.xyz)

    def test(self):
        "Test the read speed of the hdf5 code (10000 frames, 100 atoms)"
        with hdf5.HDF5TrajectoryFile(self.fn) as f:
            f.read()
