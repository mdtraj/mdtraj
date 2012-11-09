import subprocess
from distutils.core import *

# rerun swig manually. This makes sure that src/xtc.py exists. If setuptools
# runs build before build_ext, then without this line the xtc.py file will
# not get moved to the right build directory and installed
subprocess.check_call('swig -python -c++ -o src/xtc_wrap.cpp src/xtc.i', shell=True)

xtc = Extension("_xtc",
                ["src/xtc_wrap.cpp","src/xtc.cpp",
                 'xdrlib/xdrfile.c',
                 'xdrlib/trr2xtc.c',
                 'xdrlib/xdrfile_trr.c',
                 'xdrlib/xdrfile_xtc.c'],
                include_dirs=['src', 'xdrlib/include'],
                )


# NumyTypemapTests setup
setup(name = "err test",
      py_modules = ['xtc'],
      package_dir = {'' : 'src'},
      ext_modules = [xtc])

