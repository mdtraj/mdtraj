"""MDTraj: a python library for loading, saving, and manipulating molecular dynamics trajectories.

MDTraj provides an easy to use python interface for manipulating MD trajectories.
It supports the reading and writing of molecular dynamics trajectories
in a variety of formats, including full support for PDB, DCD, XTC, TRR, binpos,
AMBER NetCDF, and MDTraj HDF5. The package also provides a command line script
for converting trajectories between supported formats.
"""

DOCLINES = __doc__.split("\n")

import os
import subprocess
from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy

##########################
VERSION = "0.3.2"
ISRELEASED = False
__version__ = VERSION
##########################

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Programming Language :: C
Programming Language :: Python
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

rmsd = Extension('mdtraj._rmsd',
    sources = ['MDTraj/rmsd/src/theobald_rmsd.c','MDTraj/rmsd/_rmsd.pyx'],
    include_dirs = ["MDTraj/rmsd/include", numpy.get_include()],
    extra_compile_args = ["-std=c99","-O2", "-msse2","-msse3","-fopenmp"],
    extra_link_args = ['-lgomp', '-lm'],
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
      author='Robert McGibbon',
      author_email='rmcgibbo@gmail.com',
      description=DOCLINES[0],
      long_description="\n".join(DOCLINES[2:]),
      version=__version__,
      license='GPLv3+',
      url = "http://rmcgibbo.github.io/mdtraj",
      platforms = ["Linux", "Mac OS-X", "Unix"],
      classifiers = CLASSIFIERS.splitlines(),
      packages=['mdtraj', 'mdtraj.pdb', 'mdtraj.testing', 'mdtraj.utils',
                'mdtraj.reporters', 'mdtraj.geometry', 'mdtraj.tests'],
      package_dir={'mdtraj':'MDTraj'},
      install_requires=['numpy', 'cython', 'nose', 'nose-exclude'],
      zip_safe=False,
      scripts=['scripts/mdconvert', 'scripts/mdinspect'],
      ext_modules=[xtc, trr, dcd, binpos, rmsd],
      cmdclass = {'build_ext': build_ext},
      package_data = {'mdtraj.pdb': ['data/*'],
                      'mdtraj.testing': ["reference/*"]})

