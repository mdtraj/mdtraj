"""MDTraj: A modern, open library for the analysis of molecular dynamics trajectories

MDTraj is a python library that allows users to manipulate molecular dynamics
(MD) trajectories and perform a variety of analyses, including fast RMSD,
solvent accessible surface area, hydrogen bonding, etc. A highlight of MDTraj
is the wide variety of molecular dynamics trajectory file formats which are
supported, including RCSB pdb, GROMACS xtc, tng, and trr, CHARMM / NAMD dcd, AMBER
binpos, AMBER NetCDF, AMBER mdcrd, TINKER arc and MDTraj HDF5.
"""
from __future__ import print_function, absolute_import

import sys
from glob import glob

DOCLINES = __doc__.split("\n")

from setuptools import setup, Extension, find_packages

sys.path.insert(0, '.')
from basesetup import (write_version_py, build_ext,
                       StaticLibrary, CompilerDetection, parse_setuppy_commands)

try:
    # add an optional --disable-openmp to disable OpenMP support
    sys.argv.remove('--disable-openmp')
    disable_openmp = True
except ValueError:
    disable_openmp = False


##########################
VERSION = "1.9.9"
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
compiler.initialize()

extra_cpp_libraries = []

if sys.platform == 'win32':
    extra_cpp_libraries.append('Ws2_32')
    # For determining if a path is relative (for dtr)
    extra_cpp_libraries.append('Shlwapi')


################################################################################
# Declaration of the compiled extension modules (cython + c)
################################################################################

def format_extensions():
    compiler_args = compiler.compiler_args_warn

    xtc = Extension('mdtraj.formats.xtc',
                    sources=['mdtraj/formats/xtc/src/xdrfile.c',
                             'mdtraj/formats/xtc/src/xdr_seek.c',
                             'mdtraj/formats/xtc/src/xdrfile_xtc.c',
                             'mdtraj/formats/xtc/xtc.pyx',
                             ],
                    include_dirs=['mdtraj/formats/xtc/include/',
                                  'mdtraj/formats/xtc/'],
                    extra_compile_args=compiler_args)

    trr = Extension('mdtraj.formats.trr',
                    sources=['mdtraj/formats/xtc/src/xdrfile.c',
                             'mdtraj/formats/xtc/src/xdr_seek.c',
                             'mdtraj/formats/xtc/src/xdrfile_trr.c',
                             'mdtraj/formats/xtc/trr.pyx'],
                    include_dirs=['mdtraj/formats/xtc/include/',
                                  'mdtraj/formats/xtc/'],
                    extra_compile_args=compiler_args)

    zlib_include_dirs = []
    zlib_library_dirs = []
    if sys.platform == 'win32':
        # Conda puts the zlib headers in ./Library/... on windows
        # If you're not using conda, good luck!
        zlib_include_dirs += ["{}/Library/include".format(sys.prefix)]
        zlib_library_dirs += ["{}/Library/lib".format(sys.prefix)]
    else:
        # On linux (and mac(?)) these paths should work for a standard
        # install of python+zlib or a conda install of python+zlib
        zlib_include_dirs += ["{}/include".format(sys.prefix)]
        zlib_library_dirs += ["{}/lib".format(sys.prefix)]

    tng = Extension('mdtraj.formats.tng',
                    sources=glob('mdtraj/formats/tng/src/compression/*.c') +
                                ['mdtraj/formats/tng/src/lib/tng_io.c',
                                 'mdtraj/formats/tng/src/lib/md5.c',
                                 'mdtraj/formats/tng/tng.pyx'],
                    include_dirs=['mdtraj/formats/tng/include']
                                 + zlib_include_dirs,
                    define_macros=[('USE_ZLIB', 1)],
                    library_dirs=zlib_library_dirs,
                    libraries=['z'],
                    )

    dcd = Extension('mdtraj.formats.dcd',
                    sources=['mdtraj/formats/dcd/src/dcdplugin.c',
                             'mdtraj/formats/dcd/dcd.pyx'],
                    include_dirs=["mdtraj/formats/dcd/include/",
                                  'mdtraj/formats/dcd/'],
                    extra_compile_args=compiler_args)

    binpos = Extension('mdtraj.formats.binpos',
                       sources=['mdtraj/formats/binpos/src/binposplugin.c',
                                'mdtraj/formats/binpos/binpos.pyx'],
                       include_dirs=['mdtraj/formats/binpos/include/',
                                     'mdtraj/formats/binpos/'],
                       extra_compile_args=compiler_args)

    dtr = Extension('mdtraj.formats.dtr',
                    sources=['mdtraj/formats/dtr/src/dtrplugin.cxx',
                             'mdtraj/formats/dtr/dtr.pyx'],
                    include_dirs=['mdtraj/formats/dtr/include/',
                                  'mdtraj/formats/dtr/'],
                    define_macros=[('DESRES_READ_TIMESTEP2', 1)],
                    language='c++',
                    extra_compile_args=compiler_args,
                    libraries=extra_cpp_libraries)

    return [xtc, trr, tng, dcd, binpos, dtr]


def rmsd_extensions():
    compiler_args = (compiler.compiler_args_openmp + compiler.compiler_args_sse2 +
                     compiler.compiler_args_sse3 + compiler.compiler_args_opt +
                     compiler.compiler_args_warn)
    compiler_libraries = compiler.compiler_libraries_openmp
    define_macros = [("__NO_INTRINSICS",1)] if compiler.disable_intrinsics else None

    libtheobald = StaticLibrary(
        'mdtraj.core.lib.libtheobald',
        sources=[
            'mdtraj/rmsd/src/theobald_rmsd.cpp',
            'mdtraj/rmsd/src/center.cpp'],
        include_dirs=[
            'mdtraj/rmsd/include'],
        export_include=['mdtraj/rmsd/include/theobald_rmsd.h',
                        'mdtraj/rmsd/include/center.h'],
        language="c++",
        # don't enable OpenMP
        extra_compile_args=(compiler.compiler_args_sse2 +
                            compiler.compiler_args_sse3 +
                            compiler.compiler_args_opt))

    rmsd = Extension('mdtraj._rmsd',
                     sources=[
                         'mdtraj/rmsd/src/theobald_rmsd.cpp',
                         'mdtraj/rmsd/src/rotation.cpp',
                         'mdtraj/rmsd/src/center.cpp',
                         'mdtraj/rmsd/_rmsd.pyx'],
                     include_dirs=['mdtraj/rmsd/include'],
                     define_macros=define_macros,
                     extra_compile_args=compiler_args,
                     libraries=compiler_libraries,
                     language="c++")

    lprmsd = Extension('mdtraj._lprmsd',
                       sources=[
                           'mdtraj/rmsd/src/theobald_rmsd.cpp',
                           'mdtraj/rmsd/src/rotation.cpp',
                           'mdtraj/rmsd/src/center.cpp',
                           'mdtraj/rmsd/src/fancy_index.cpp',
                           'mdtraj/rmsd/src/Munkres.cpp',
                           'mdtraj/rmsd/src/euclidean_permutation.cpp',
                           'mdtraj/rmsd/_lprmsd.pyx'],
                       language='c++',
                       define_macros=define_macros,
                       include_dirs=['mdtraj/rmsd/include'],
                       extra_compile_args=compiler_args,
                       libraries=compiler_libraries + extra_cpp_libraries)
    return rmsd, lprmsd, libtheobald


def geometry_extensions():
    compiler.initialize()
    compiler_args = (
        compiler.compiler_args_openmp +
        compiler.compiler_args_sse2 + compiler.compiler_args_sse3 +
        compiler.compiler_args_opt + compiler.compiler_args_warn)
    define_macros = [("__NO_INTRINSICS",1)] if compiler.disable_intrinsics else None
    compiler_libraries = compiler.compiler_libraries_openmp + extra_cpp_libraries

    return [
        Extension('mdtraj.geometry._geometry',
            sources=['mdtraj/geometry/src/sasa.cpp',
                     'mdtraj/geometry/src/dssp.cpp',
                     'mdtraj/geometry/src/geometry.cpp',
                     'mdtraj/geometry/src/_geometry.pyx',],
            include_dirs=['mdtraj/geometry/include',
                          'mdtraj/geometry/src/kernels'],
            depends=['mdtraj/geometry/src/kernels/anglekernels.h',
                     'mdtraj/geometry/src/kernels/dihedralkernels.h',
                     'mdtraj/geometry/src/kernels/distancekernels.h'],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            libraries=compiler_libraries,
            language='c++'),
        Extension('mdtraj.geometry.drid',
            sources=["mdtraj/geometry/drid.pyx",
                     "mdtraj/geometry/src/dridkernels.cpp",
                     "mdtraj/geometry/src/moments.cpp"],
            include_dirs=["mdtraj/geometry/include"],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            libraries=compiler_libraries,
            language='c++'),
        Extension('mdtraj.geometry.neighbors',
            sources=["mdtraj/geometry/neighbors.pyx",
                     "mdtraj/geometry/src/neighbors.cpp"],
            include_dirs=["mdtraj/geometry/include",],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            libraries=compiler_libraries,
            language='c++'),
        Extension('mdtraj.geometry.neighborlist',
            sources=["mdtraj/geometry/neighborlist.pyx",
                     "mdtraj/geometry/src/neighborlist.cpp"],
            include_dirs=["mdtraj/geometry/include",],
            define_macros=define_macros,
            extra_compile_args=compiler_args,
            libraries=compiler_libraries,
            language='c++'),
        ]


write_version_py(VERSION, ISRELEASED, 'mdtraj/version.py')

metadata = \
    dict(name='mdtraj',
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
      install_requires=['numpy>=1.6',
                        'scipy',
                        'astunparse',
                        'pyparsing',
                        ],
      package_data={'mdtraj.formats.pdb': ['data/*'], },
      zip_safe=False,
      entry_points={'console_scripts':
          ['mdconvert = mdtraj.scripts.mdconvert:entry_point',
           'mdinspect = mdtraj.scripts.mdinspect:entry_point']},
)


if __name__ == '__main__':
    # Don't use numpy if we are just - non-build actions are required to succeed
    # without NumPy for example when pip is used to install Scipy when
    # NumPy is not yet present in the system.
    run_build = parse_setuppy_commands()
    if run_build:
        extensions = format_extensions()
        extensions.extend(rmsd_extensions())
        extensions.extend(geometry_extensions())

        # most extensions use numpy, add headers for it.
        try:
            import Cython as _c
            from Cython.Build import cythonize
            if _c.__version__ < '0.29':
                raise ImportError("Too old")
        except ImportError as e:
            print('mdtrajs setup depends on Cython (>=0.29). Install it prior invoking setup.py')
            print(e)
            sys.exit(1)
        try:
            import numpy as np
        except ImportError:
            print('mdtrajs setup depends on NumPy. Install it prior invoking setup.py')
            sys.exit(1)

        for e in extensions:
            e.include_dirs.append(np.get_include())
        metadata['ext_modules'] = cythonize(
            extensions,
            language_level=sys.version_info[0],
            force=True,
        )

    setup(**metadata)
