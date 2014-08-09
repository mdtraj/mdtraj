"""mdtraj: Read, write and analyze MD trajectories with only a few lines of Python code.

mdtraj is a python library that allows users to manipulate molecular dynamics
(MD) trajectories and perform a variety of analyses, including fast RMSD,
solvent accessible surface area, hydrogen bonding, etc. A highlight of mdtraj
is the wide variety of molecular dynamics trajectory file formats which are
supported, including RCSB pdb, GROMACS xtc and trr, CHARMM / NAMD dcd, AMBER
binpos, AMBER NetCDF, AMBER mdcrd, TINKER arc and mdtraj HDF5.
"""

from __future__ import print_function
DOCLINES = __doc__.split("\n")

import os
import sys
import shutil
import tempfile
import subprocess
from distutils.ccompiler import new_compiler
try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

try:
    import numpy
except ImportError:
    print('Building and running mdtraj requires numpy', file=sys.stderr)
    sys.exit(1)

try:
    import Cython
    if Cython.__version__ < '0.19':
        raise ImportError
    from Cython.Build import cythonize
except ImportError:
    print('Building from source requires cython >= 0.19', file=sys.stderr)
    exit(1)

try:
    # add an optional --disable-openmp to disable OpenMP support
    sys.argv.remove('--disable-openmp')
    disable_openmp = True
except ValueError:
    disable_openmp = False

try:
    # add an optional command line flag --no-install-deps to setup.py
    # to turn off setuptools automatic downloading of dependencies
    sys.argv.remove('--no-install-deps')
    no_install_deps = True
except ValueError:
    no_install_deps = False


setup_kwargs = {}
if 'setuptools' in sys.modules:
    setup_kwargs['zip_safe'] = False
    setup_kwargs['entry_points'] = {'console_scripts':
              ['mdconvert = mdtraj.scripts.mdconvert:entry_point',
               'mdinspect = mdtraj.scripts.mdinspect:entry_point']}

    if sys.version_info[0] == 2:
        # required to fix cythoninze() for old versions of setuptools on
        # python 2
        m = sys.modules['setuptools.extension']
        m.Extension.__dict__ = m._Extension.__dict__
else:
    setup_kwargs['scripts'] = ['scripts/mdconvert.py', 'scripts/mdinspect.py']


##########################
VERSION = "0.9.1"
ISRELEASED = False
__version__ = VERSION
##########################


CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

def find_packages():
    """Find all of mdtraj's python packages.
    Adapted from IPython's setupbase.py. Copyright IPython
    contributors, licensed under the BSD license.
    """
    packages = ['mdtraj.scripts']
    for dir,subdirs,files in os.walk('mdtraj'):
        package = dir.replace(os.path.sep, '.')
        if '__init__.py' not in files:
            # not a package
            continue
        packages.append(package.replace('mdtraj', 'mdtraj'))
    return packages


################################################################################
# Writing version control information to the module
################################################################################

def git_version():
    # Return the git revision as a string
    # copied from numpy setup.py
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
        out = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = 'Unknown'

    return GIT_REVISION


def write_version_py(filename='mdtraj/version.py'):
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
        GIT_REVISION = 'Unknown'

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

################################################################################
# Detection of compiler capabilities
################################################################################

class CompilerDetection(object):
    def __init__(self, disable_openmp):
        self.msvc = new_compiler().compiler_type == 'msvc'

        if disable_openmp:
            self.openmp_enabled = False
        else:
            self.openmp_enabled, openmp_needs_gomp = self._detect_openmp()
        self.sse3_enabled = self._detect_sse3() if not self.msvc else True
        self.sse41_enabled = self._detect_sse41() if not self.msvc else True

        self.compiler_args_sse2  = ['-msse2'] if not self.msvc else ['/arch:SSE2']
        self.compiler_args_sse3  = ['-mssse3'] if (self.sse3_enabled and not self.msvc) else []

        self.compiler_args_sse41, self.define_macros_sse41 = [], []
        if self.sse41_enabled:
            self.define_macros_sse41 = [('__SSE4__', 1), ('__SSE4_1__', 1)]
            if not self.msvc:
                self.compiler_args_sse41 = ['-msse4']

        if self.openmp_enabled:
            self.compiler_libraries_openmp = ['gomp'] if openmp_needs_gomp else []

            if self.msvc:
                self.compiler_args_openmp = ['/openmp']
            else:
                self.compiler_args_openmp = ['-fopenmp']
        else:
            self.compiler_libraries_openmp = []
            self.compiler_args_openmp = []

        if self.msvc:
            self.compiler_args_opt = ['/O2']
        else:
            self.compiler_args_opt = ['-O3', '-funroll-loops']

    def hasfunction(self, cc, funcname, include=None, extra_postargs=None):
        # From http://stackoverflow.com/questions/
        #            7018879/disabling-output-when-compiling-with-distutils
        tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
        devnull = oldstderr = None
        try:
            try:
                fname = os.path.join(tmpdir, 'funcname.c')
                f = open(fname, 'w')
                if include is not None:
                    f.write('#include %s\n' % include)
                f.write('int main(void) {\n')
                f.write('    %s;\n' % funcname)
                f.write('}\n')
                f.close()
                devnull = open(os.devnull, 'w')
                oldstderr = os.dup(sys.stderr.fileno())
                os.dup2(devnull.fileno(), sys.stderr.fileno())
                objects = cc.compile([fname], output_dir=tmpdir,
                                     extra_postargs=extra_postargs)
                cc.link_executable(objects, os.path.join(tmpdir, 'a.out'))
            except Exception as e:
                return False
            return True
        finally:
            if oldstderr is not None:
                os.dup2(oldstderr, sys.stderr.fileno())
            if devnull is not None:
                devnull.close()
            shutil.rmtree(tmpdir)

    def _print_support_start(self, feature):
        print('Attempting to autodetect {0} support...'.format(feature), end=' ')

    def _print_support_end(self, feature, status):
        if status is True:
            print('Compiler supports {0}'.format(feature))
        else:
            print('Did not detect {0} support'.format(feature))

    def _detect_openmp(self):
        self._print_support_start('OpenMP')
        compiler = new_compiler()
        hasopenmp = self.hasfunction(compiler, 'omp_get_num_threads()', extra_postargs=['-fopenmp', '/openmp'])
        needs_gomp = hasopenmp
        if not hasopenmp:
            compiler.add_library('gomp')
            hasopenmp = self.hasfunction(compiler, 'omp_get_num_threads()')
            needs_gomp = hasopenmp
        self._print_support_end('OpenMP', hasopenmp)
        return hasopenmp, needs_gomp

    def _detect_sse3(self):
        "Does this compiler support SSE3 intrinsics?"
        compiler = new_compiler()
        self._print_support_start('SSE3')
        result = self.hasfunction(compiler, '__m128 v; _mm_hadd_ps(v,v)',
                           include='<pmmintrin.h>',
                           extra_postargs=['-msse3'])
        self._print_support_end('SSE3', result)
        return result

    def _detect_sse41(self):
        "Does this compiler support SSE4.1 intrinsics?"
        compiler = new_compiler()
        self._print_support_start('SSE4.1')
        result = self.hasfunction(compiler, '__m128 v; _mm_round_ps(v,0x00)',
                           include='<smmintrin.h>',
                           extra_postargs=['-msse4'])
        self._print_support_end('SSE4.1', result)
        return result

# Global info about compiler
compiler = CompilerDetection(disable_openmp)

################################################################################
# Declaration of the compiled extension modules (cython + c)
################################################################################


xtc = Extension('mdtraj.formats.xtc',
                sources=['mdtraj/formats/xtc/src/xdrfile.c',
                         'mdtraj/formats/xtc/src/xdrfile_xtc.c',
                         'mdtraj/formats/xtc/xtc.pyx'],
                include_dirs=['mdtraj/formats/xtc/include/',
                              'mdtraj/formats/xtc/', numpy.get_include()])

trr = Extension('mdtraj.formats.trr',
                sources=['mdtraj/formats/xtc/src/xdrfile.c',
                         'mdtraj/formats/xtc/src/xdrfile_trr.c',
                         'mdtraj/formats/xtc/trr.pyx'],
                include_dirs=['mdtraj/formats/xtc/include/',
                              'mdtraj/formats/xtc/', numpy.get_include()])

dcd = Extension('mdtraj.formats.dcd',
                sources=['mdtraj/formats/dcd/src/dcdplugin.c',
                         'mdtraj/formats/dcd/dcd.pyx'],
                include_dirs=["mdtraj/formats/dcd/include/",
                              'mdtraj/formats/dcd/', numpy.get_include()])

binpos = Extension('mdtraj.formats.binpos',
                   sources=['mdtraj/formats/binpos/src/binposplugin.c',
                            'mdtraj/formats/binpos/binpos.pyx'],
                   include_dirs=['mdtraj/formats/binpos/include/',
                                 'mdtraj/formats/binpos/', numpy.get_include()])


def rmsd_extensions():
    compiler_args = (compiler.compiler_args_openmp + compiler.compiler_args_sse2 +
                     compiler.compiler_args_sse3 + compiler.compiler_args_opt)
    compiler_libraries = compiler.compiler_libraries_openmp
    rmsd = Extension('mdtraj._rmsd',
                     sources=[
                         'mdtraj/rmsd/src/theobald_rmsd.c',
                         'mdtraj/rmsd/src/rotation.c',
                         'mdtraj/rmsd/src/center.c',
                         'mdtraj/rmsd/_rmsd.pyx'],
                     include_dirs=[
                         'mdtraj/rmsd/include', numpy.get_include()],
                     extra_compile_args=compiler_args,
                     libraries=compiler_libraries)

    lprmsd = Extension('mdtraj._lprmsd',
                       sources=[
                           'mdtraj/rmsd/src/theobald_rmsd.c',
                           'mdtraj/rmsd/src/rotation.c',
                           'mdtraj/rmsd/src/center.c',
                           'mdtraj/rmsd/src/fancy_index.cpp',
                           'mdtraj/rmsd/src/Munkres.cpp',
                           'mdtraj/rmsd/src/euclidean_permutation.cpp',
                           'mdtraj/rmsd/_lprmsd.pyx'],
                       language='c++',
                       include_dirs=[
                           'mdtraj/rmsd/include', numpy.get_include()],
                       extra_compile_args=compiler_args,
                       libraries=compiler_libraries)
    return rmsd, lprmsd



def geometry_extensions():
    compiler_args = (compiler.compiler_args_sse2 + compiler.compiler_args_sse3 + compiler.compiler_args_sse41 +
                     compiler.compiler_args_opt + compiler.compiler_args_sse41)
    define_macros = compiler.define_macros_sse41

    return [
        Extension('mdtraj.geometry._geometry',
            sources=['mdtraj/geometry/src/geometry.c',
                     'mdtraj/geometry/src/sasa.c',
                     'mdtraj/geometry/src/_geometry.pyx'],
            include_dirs=['mdtraj/geometry/include',
                          'mdtraj/geometry/src/kernels',
                          numpy.get_include()],
            define_macros=define_macros,
            extra_compile_args=compiler_args),
        Extension('mdtraj.geometry.drid',
            sources=["mdtraj/geometry/drid.pyx",
                     "mdtraj/geometry/src/dridkernels.c",
                     "mdtraj/geometry/src/cephes/cbrt.c",
                     "mdtraj/geometry/src/cephes/isnan.c",
                     "mdtraj/geometry/src/moments.c"],
            include_dirs=["mdtraj/geometry/include",
                          "mdtraj/geometry/include/cephes",
                          numpy.get_include()],
            define_macros=define_macros,
            extra_compile_args=compiler_args)
        ]

extensions = [xtc, trr, dcd, binpos]
extensions.extend(rmsd_extensions())
extensions.extend(geometry_extensions())

write_version_py()
setup(name='mdtraj',
      author='Robert McGibbon',
      author_email='rmcgibbo@gmail.com',
      description=DOCLINES[0],
      long_description="\n".join(DOCLINES[2:]),
      version=__version__,
      license='LGPLv2.1+',
      url='http://mdtraj.org',
      download_url = "https://github.com/rmcgibbo/mdtraj/releases/latest",
      platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
      classifiers=CLASSIFIERS.splitlines(),
      packages=find_packages(),
      package_dir={'mdtraj': 'mdtraj', 'mdtraj.scripts': 'scripts'},
      ext_modules=cythonize(extensions),
      package_data={'mdtraj.formats.pdb': ['data/*'],
                    'mdtraj.testing': ['reference/*']},
      **setup_kwargs)
