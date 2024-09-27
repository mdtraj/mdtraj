import numpy as np
from packaging.version import Version

# This is a debug mode used to catch changes introduced by NEP50.
# This should be deleted after everything is fixed in mdtraj.
# See https://numpy.org/neps/nep-0050-scalar-promotion.html
if Version(np.__version__) >= Version('2.0.0a0'):
    np._set_promotion_state("weak_and_warn")
