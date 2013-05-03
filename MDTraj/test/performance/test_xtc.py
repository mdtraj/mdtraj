import os
import tempfile
from noseperf.testcases import PerformanceTest

from mdtraj import xtc
import numpy as np

######################################################
# Base class: handles setting up a temp file and array
######################################################

class WithTemp(PerformanceTest):
    n_frames = 10000
    n_atoms = 100
    
    def setUp(self):
        self.xyz = np.random.randn(self.n_frames, self.n_atoms, 3)
        self.fn = tempfile.mkstemp()[1]

    def tearDown(self):
        os.unlink(self.fn)

########################################
# Tests
########################################

class TestWrite(WithTemp):
    def test(self):
        "Test the write speed of the XTC code"
        xtc.write(self.fn, xyz=self.xyz, force_overwrite=True)

class TestRead(WithTemp):
    def setUp(self):
        super(TestRead, self).setUp()
        xtc.write(self.fn, xyz=self.xyz, force_overwrite=True)

    def test(self):
        "Test the read speed of the XTC code"
        xtc.read(self.fn)
