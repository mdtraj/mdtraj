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

import os
import subprocess
from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy

##########################
VERSION = "0.2"
ISRELEASED = False
__version__ = VERSION
##########################


IRMSD = Extension('mdtraj.IRMSD',
    sources = ['MDTraj/IRMSD/theobald_rmsd.c','MDTraj/IRMSD/rmsd_wrap.pyx'],
    include_dirs = ["MDTraj/IRMSD/", numpy.get_include()],
    extra_compile_args = ["-std=c99","-O2", "-msse2","-msse3","-fopenmp"],
    extra_link_args = ['-lgomp'],
    )

xtc = Extension('mdtraj.xtc',
    sources = ['MDTraj/xtc/src/xdrfile.c', 'MDTraj/xtc/src/xdrfile_xtc.c',
               'MDTraj/xtc/xtc.pyx'],
    include_dirs = ["MDTraj/xtc/include/", 'MDTraj/xtc/', numpy.get_include()])

trr = Extension('mdtraj.trr',
    sources = ['MDTraj/xtc/src/xdrfile.c', 'MDTraj/xtc/src/xdrfile_trr.c', 
               'MDTraj/xtc/trr.pyx'],
    include_dirs = ["MDTraj/xtc/include/", 'MDTraj/xtc/', numpy.get_include()])

dcd = Extension('mdtraj.dcd',
    sources = ["MDTraj/dcd/src/dcdplugin.c", "MDTraj/dcd/dcd.pyx"],
    libraries=['m'],
    include_dirs = ["MDTraj/dcd/include/", 'MDTraj/dcd/', numpy.get_include()])

binpos = Extension('mdtraj.binpos',
    sources = ['MDTraj/binpos/src/binposplugin.c', 'MDTraj/binpos/binpos.pyx'],
    include_dirs = ["MDTraj/binpos/include/", 'MDTraj/binpos/', numpy.get_include()])


# Return the git revision as a string
# copied from numpy setup.py
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

def write_version_py(filename='MDTraj/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM MDTRAJ SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of numpy.version messes up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


write_version_py()
setup(name='mdtraj',
      version=__version__,
      packages=['mdtraj', 'mdtraj.pdb', 'mdtraj.testing', 'mdtraj.utils',
                'mdtraj.reporters', 'mdtraj.geometry'],
      package_dir={'mdtraj':'MDTraj'},
      install_requires=['numpy', 'cython', 'nose', 'nose-exclude'],
      zip_safe=False,
      scripts=['scripts/mdconvert'],
      ext_modules=[xtc, trr, dcd, binpos, IRMSD],
      cmdclass = {'build_ext': build_ext},
      package_data = {'mdtraj.pdb': ['data/*'],
                      'mdtraj.testing': ["reference/*"]})

