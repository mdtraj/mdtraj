"""MDTraj: A modern, open library for the analysis of molecular dynamics trajectories

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


if sys.version_info[0] == 2:
    # required to fix cythoninze() for old versions of setuptools on
    # python 2
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__


##########################
VERSION = "1.4.2"
ISRELEASED = True
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

dtr = Extension('mdtraj.formats.dtr',
                   sources=['mdtraj/formats/dtr/src/dtrplugin.cxx',
                            'mdtraj/formats/dtr/dtr.pyx'],
                   include_dirs=['mdtraj/formats/dtr/include/',
                                 'mdtraj/formats/dtr/', numpy.get_include()],
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
            'mdtraj/rmsd/src/theobald_rmsd.c',
            'mdtraj/rmsd/src/center.c'],
        include_dirs=[
            'mdtraj/rmsd/include'],
        export_include=['mdtraj/rmsd/include/theobald_rmsd.h',
                        'mdtraj/rmsd/include/center.h'],
        # don't enable OpenMP
        extra_compile_args=(compiler.compiler_args_sse2 +
                            compiler.compiler_args_sse3 +
                            compiler.compiler_args_opt))

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
                       libraries=compiler_libraries + extra_cpp_libraries)
    return rmsd, lprmsd, libtheobald



def geometry_extensions():
    compiler_args = (compiler.compiler_args_sse2 + compiler.compiler_args_sse3 +
                     compiler.compiler_args_opt)
    define_macros = None

    return [
        Extension('mdtraj.geometry._geometry',
            sources=['mdtraj/geometry/src/geometry.c',
                     'mdtraj/geometry/src/sasa.c',
                     'mdtraj/geometry/src/dssp.cpp',
                     'mdtraj/geometry/src/_geometry.pyx'],
            include_dirs=['mdtraj/geometry/include',
                          'mdtraj/geometry/src/kernels',
                          numpy.get_include()],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            libraries=extra_cpp_libraries,
            language='c++'),
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
            extra_compile_args=compiler_args),
        Extension('mdtraj.geometry.neighbors',
            sources=["mdtraj/geometry/neighbors.pyx",
                     "mdtraj/geometry/src/neighbors.cpp"],
            include_dirs=["mdtraj/geometry/include",],
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
      zip_safe=False,
      entry_points={'console_scripts':
          ['mdconvert = mdtraj.scripts.mdconvert:entry_point',
           'mdinspect = mdtraj.scripts.mdinspect:entry_point']},
)
