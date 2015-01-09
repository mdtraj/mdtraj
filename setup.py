"""MDTraj: Read, write and analyze MD trajectories with only a few lines of Python code.

MDTraj is a python library that allows users to manipulate molecular dynamics
(MD) trajectories and perform a variety of analyses, including fast RMSD,
solvent accessible surface area, hydrogen bonding, etc. A highlight of MDTraj
is the wide variety of molecular dynamics trajectory file formats which are
supported, including RCSB pdb, GROMACS xtc and trr, CHARMM / NAMD dcd, AMBER
binpos, AMBER NetCDF, AMBER mdcrd, TINKER arc and MDTraj HDF5.
"""

from __future__ import print_function, absolute_import
DOCLINES = __doc__.split("\n")

import sys
from setuptools import setup, Extension
sys.path.insert(0, '.')
from basesetup import (find_packages, write_version_py, build_ext,
                       StaticLibrary, CompilerDetection)

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
VERSION = "1.3.0.dev0"
ISRELEASED = False
__version__ = VERSION
##########################


CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
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

# Global info about compiler
compiler = CompilerDetection(disable_openmp)

extra_cpp_libraries = []
if sys.platform == 'darwin':
    extra_cpp_libraries.append('stdc++')
if sys.platform == 'win32':
    extra_cpp_libraries.append('Ws2_32')
    # For determining if a path is relative (for dtr)
    extra_cpp_libraries.append('Shlwapi')


################################################################################
# Declaration of the compiled extension modules (cython + c)
################################################################################


xtc = Extension('mdtraj.formats.xtc',
                sources=['MDTraj/formats/xtc/src/xdrfile.c',
                         'MDTraj/formats/xtc/src/xdrfile_xtc.c',
                         'MDTraj/formats/xtc/xtc.pyx'],
                include_dirs=['MDTraj/formats/xtc/include/',
                              'MDTraj/formats/xtc/', numpy.get_include()])

trr = Extension('mdtraj.formats.trr',
                sources=['MDTraj/formats/xtc/src/xdrfile.c',
                         'MDTraj/formats/xtc/src/xdrfile_trr.c',
                         'MDTraj/formats/xtc/trr.pyx'],
                include_dirs=['MDTraj/formats/xtc/include/',
                              'MDTraj/formats/xtc/', numpy.get_include()])

dcd = Extension('mdtraj.formats.dcd',
                sources=['MDTraj/formats/dcd/src/dcdplugin.c',
                         'MDTraj/formats/dcd/dcd.pyx'],
                include_dirs=["MDTraj/formats/dcd/include/",
                              'MDTraj/formats/dcd/', numpy.get_include()])

binpos = Extension('mdtraj.formats.binpos',
                   sources=['MDTraj/formats/binpos/src/binposplugin.c',
                            'MDTraj/formats/binpos/binpos.pyx'],
                   include_dirs=['MDTraj/formats/binpos/include/',
                                 'MDTraj/formats/binpos/', numpy.get_include()])

dtr = Extension('mdtraj.formats.dtr',
                   sources=['MDTraj/formats/dtr/src/dtrplugin.cxx',
                            'MDTraj/formats/dtr/dtr.pyx'],
                   include_dirs=['MDTraj/formats/dtr/include/',
                                 'MDTraj/formats/dtr/', numpy.get_include()],
                   define_macros = [('DESRES_READ_TIMESTEP2', 1)],
                   language='c++',
                   libraries=extra_cpp_libraries)


def rmsd_extensions():
    compiler_args = (compiler.compiler_args_openmp + compiler.compiler_args_sse2 +
                     compiler.compiler_args_sse3 + compiler.compiler_args_opt)
    compiler_libraries = compiler.compiler_libraries_openmp

    libtheobald = StaticLibrary(
        'mdtraj.core.lib.libtheobald',
        sources=[
            'MDTraj/rmsd/src/theobald_rmsd.c',
            'MDTraj/rmsd/src/center.c'],
        include_dirs=[
            'MDTraj/rmsd/include'],
        export_include=['MDTraj/rmsd/include/theobald_rmsd.h',
                        'MDTraj/rmsd/include/center.h'],
        # don't enable OpenMP
        extra_compile_args=(compiler.compiler_args_sse2 +
                            compiler.compiler_args_sse3 +
                            compiler.compiler_args_opt))

    rmsd = Extension('mdtraj._rmsd',
                     sources=[
                         'MDTraj/rmsd/src/theobald_rmsd.c',
                         'MDTraj/rmsd/src/rotation.c',
                         'MDTraj/rmsd/src/center.c',
                         'MDTraj/rmsd/_rmsd.pyx'],
                     include_dirs=[
                         'MDTraj/rmsd/include', numpy.get_include()],
                     extra_compile_args=compiler_args,
                     libraries=compiler_libraries)

    lprmsd = Extension('mdtraj._lprmsd',
                       sources=[
                           'MDTraj/rmsd/src/theobald_rmsd.c',
                           'MDTraj/rmsd/src/rotation.c',
                           'MDTraj/rmsd/src/center.c',
                           'MDTraj/rmsd/src/fancy_index.cpp',
                           'MDTraj/rmsd/src/Munkres.cpp',
                           'MDTraj/rmsd/src/euclidean_permutation.cpp',
                           'MDTraj/rmsd/_lprmsd.pyx'],
                       language='c++',
                       include_dirs=[
                           'MDTraj/rmsd/include', numpy.get_include()],
                       extra_compile_args=compiler_args,
                       libraries=compiler_libraries + extra_cpp_libraries)
    return rmsd, lprmsd, libtheobald



def geometry_extensions():
    compiler_args = (compiler.compiler_args_sse2 + compiler.compiler_args_sse3 +
                     compiler.compiler_args_opt)
    define_macros = None

    return [
        Extension('mdtraj.geometry._geometry',
            sources=['MDTraj/geometry/src/geometry.c',
                     'MDTraj/geometry/src/sasa.c',
                     'MDTraj/geometry/src/dssp.cpp',
                     'MDTraj/geometry/src/_geometry.pyx'],
            include_dirs=['MDTraj/geometry/include',
                          'MDTraj/geometry/src/kernels',
                          numpy.get_include()],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            libraries=extra_cpp_libraries,
            language='c++'),
        Extension('mdtraj.geometry.drid',
            sources=["MDTraj/geometry/drid.pyx",
                     "MDTraj/geometry/src/dridkernels.c",
                     "MDTraj/geometry/src/cephes/cbrt.c",
                     "MDTraj/geometry/src/cephes/isnan.c",
                     "MDTraj/geometry/src/moments.c"],
            include_dirs=["MDTraj/geometry/include",
                          "MDTraj/geometry/include/cephes",
                          numpy.get_include()],
            define_macros=define_macros,
            extra_compile_args=compiler_args),
        Extension('mdtraj.geometry.neighbors',
            sources=["MDTraj/geometry/neighbors.pyx",
                     "MDTraj/geometry/src/neighbors.cpp"],
            include_dirs=["MDTraj/geometry/include",],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            language='c++'),
        ]

extensions = [xtc, trr, dcd, binpos, dtr]
extensions.extend(rmsd_extensions())
extensions.extend(geometry_extensions())

write_version_py(VERSION, ISRELEASED)
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
      package_dir={'mdtraj': 'MDTraj', 'mdtraj.scripts': 'scripts'},
      cmdclass={'build_ext': build_ext},
      ext_modules=cythonize(extensions),
      package_data={'mdtraj.formats.pdb': ['data/*'],
                    'mdtraj.testing': ['reference/*',
                                       'reference/ala_dipeptide_trj/*',
                                       'reference/ala_dipeptide_trj/not_hashed/*',
                                       'reference/frame0.dtr/*',
                                       'reference/frame0.dtr/not_hashed/*',],
                    'mdtraj.html': ['static/*']},
      exclude_package_data={'mdtraj.testing': ['reference/ala_dipeptide_trj',
                                               'reference/ala_dipeptide_trj/not_hashed',
                                               'reference/frame0.dtr',
                                               'reference/frame0.dtr/not_hashed',]},
      **setup_kwargs)
