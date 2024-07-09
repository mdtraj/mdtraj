import json
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from distutils.ccompiler import new_compiler
from distutils.dep_util import newer_group
from distutils.errors import DistutilsExecError, DistutilsSetupError
from distutils.sysconfig import customize_compiler, get_config_vars

from setuptools import Extension
from setuptools.command.build_ext import build_ext as _build_ext

################################################################################
# Detection of compiler capabilities
################################################################################


class CompilerDetection:
    # Necessary for OSX. See https://github.com/mdtraj/mdtraj/issues/576
    # The problem is that distutils.sysconfig.customize_compiler()
    # is necessary to properly invoke the correct compiler for this class
    # (otherwise the CC env variable isn't respected). Unfortunately,
    # distutils.sysconfig.customize_compiler() DIES on OSX unless some
    # appropriate initialization routines have been called. This line
    # has a side effect of calling those initialzation routes, and is therefor
    # necessary for OSX, even though we don't use the result.
    _DONT_REMOVE_ME = get_config_vars()

    def __init__(self, disable_openmp):
        self.disable_openmp = disable_openmp
        self._is_initialized = False

    def initialize(self):
        if self._is_initialized:
            return

        cc = new_compiler()
        customize_compiler(cc)

        self.msvc = cc.compiler_type == "msvc"
        self._print_compiler_version(cc)

        if self.disable_openmp:
            self.openmp_enabled = False
        else:
            self.openmp_enabled, openmp_needs_gomp = self._detect_openmp()
        self.sse2_enabled = self._detect_sse2() if not self.msvc else True
        self.sse3_enabled = self._detect_sse3() if not self.msvc else True
        self.sse41_enabled = self._detect_sse41() if not self.msvc else True
        self.neon_enabled = self._detect_neon() if not self.msvc else False

        if self.msvc:
            self.compiler_args_sse2 = ["/arch:SSE2"]
        elif self.sse2_enabled:
            self.compiler_args_sse2 = ["-msse2"]
        else:
            self.compiler_args_sse2 = []
        self.compiler_args_sse3 = ["-mssse3"] if (self.sse3_enabled and not self.msvc) else []
        self.compiler_args_neon = []
        self.compiler_args_warn = (
            ["-Wno-unused-function", "-Wno-unreachable-code", "-Wno-sign-compare"] if not self.msvc else []
        )

        if self.neon_enabled:
            self.compiler_args_sse2 = []
            self.compiler_args_sse3 = []

        self.compiler_args_sse41, self.define_macros_sse41 = [], []
        if self.sse41_enabled:
            self.define_macros_sse41 = [("__SSE4__", 1), ("__SSE4_1__", 1)]
            if not self.msvc:
                self.compiler_args_sse41 = ["-msse4"]

        if os.environ.get("MDTRAJ_BUILD_DISABLE_INTRINSICS") == "1":
            self.disable_intrinsics = True
            print("Env var MDTRAJ_BUILD_DISABLE_INTRINSICS set")
        else:
            self.disable_intrinsics = False

        if self.openmp_enabled:
            self.compiler_libraries_openmp = []

            if self.msvc:
                self.compiler_args_openmp = ["/openmp"]
            else:
                self.compiler_args_openmp = ["-fopenmp"]
                if openmp_needs_gomp:
                    self.compiler_libraries_openmp = ["gomp"]
        else:
            self.compiler_libraries_openmp = []
            self.compiler_args_openmp = []

        if self.msvc:
            self.compiler_args_opt = ["/O2"]
        else:
            self.compiler_args_opt = ["-O3", "-funroll-loops", "--std=c++11"]
        print()
        self._is_initialized = True

    def _print_compiler_version(self, cc):
        print("C compiler:")
        try:
            if self.msvc:
                if not cc.initialized:
                    cc.initialize()
                cc.spawn([cc.cc])
            else:
                cc.spawn([cc.compiler[0]] + ["-v"])
        except DistutilsExecError:
            pass

    def hasfunction(self, funcname, include=None, libraries=None, extra_postargs=None):
        # running in a separate subshell lets us prevent unwanted stdout/stderr
        part1 = f"""
from __future__ import print_function
import os
import json
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler, get_config_vars

FUNCNAME = json.loads('{json.dumps(funcname)}')
INCLUDE = json.loads('{json.dumps(include)}')
LIBRARIES = json.loads('{json.dumps(libraries or [])}')
EXTRA_POSTARGS = json.loads('{json.dumps(extra_postargs)}')
        """

        part2 = """
get_config_vars()  # DON'T REMOVE ME
cc = new_compiler()
customize_compiler(cc)
for library in LIBRARIES:
    cc.add_library(library)

status = 0
try:
    with open('func.c', 'w') as f:
        if INCLUDE is not None:
            f.write('#include %s\\n' % INCLUDE)
        f.write('int main(void) {\\n')
        f.write('    %s;\\n' % FUNCNAME)
        f.write('}\\n')
    objects = cc.compile(['func.c'], output_dir='.',
                         extra_postargs=EXTRA_POSTARGS)
    cc.link_executable(objects, 'a.out')
except Exception as e:
    status = 1
exit(status)
        """
        tmpdir = tempfile.mkdtemp(prefix="hasfunction-")
        try:
            curdir = os.path.abspath(os.curdir)
            os.chdir(tmpdir)
            with open("script.py", "w") as f:
                f.write(part1 + part2)
            proc = subprocess.Popen(
                [sys.executable, "script.py"],
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )
            proc.communicate()
            status = proc.wait()
        finally:
            os.chdir(curdir)
            shutil.rmtree(tmpdir)

        return status == 0

    def _print_support_start(self, feature):
        print(f"Attempting to autodetect {feature:6} support...", end=" ")

    def _print_support_end(self, feature, status):
        if status is True:
            print(f"Compiler supports {feature}")
        else:
            print(f"Did not detect {feature} support")

    def _detect_openmp(self):
        self._print_support_start("OpenMP")
        extra_postargs = ["/openmp"] if self.msvc else ["-fopenmp"]
        args = dict(extra_postargs=extra_postargs, include="<omp.h>")
        hasopenmp = self.hasfunction("omp_get_num_threads()", **args)
        needs_gomp = False
        if not hasopenmp:
            hasopenmp = self.hasfunction(
                "omp_get_num_threads()",
                libraries=["gomp"],
                **args,
            )
            needs_gomp = hasopenmp
        self._print_support_end("OpenMP", hasopenmp)
        return hasopenmp, needs_gomp

    def _detect_sse2(self):
        "Does this compiler support SSE2 intrinsics?"
        self._print_support_start("SSE2")
        result = self.hasfunction(
            "__m128d v; _mm_add_pd(v,v)",
            include="<immintrin.h>",  # Possibly emmintrin.h
            extra_postargs=["-msse2"],
        )
        self._print_support_end("SSE2", result)
        return result

    def _detect_sse3(self):
        "Does this compiler support SSE3 intrinsics?"
        self._print_support_start("SSE3")
        result = self.hasfunction(
            "__m128 v; _mm_hadd_ps(v,v)",
            include="<pmmintrin.h>",
            extra_postargs=["-msse3"],
        )
        self._print_support_end("SSE3", result)
        return result

    def _detect_sse41(self):
        "Does this compiler support SSE4.1 intrinsics?"
        self._print_support_start("SSE4.1")
        result = self.hasfunction(
            "__m128 v; _mm_round_ps(v,0x00)",
            include="<smmintrin.h>",
            extra_postargs=["-msse4"],
        )
        self._print_support_end("SSE4.1", result)
        return result

    def _detect_neon(self):
        """Does this compiler support NEON intrinsics (ARM64)"""
        self._print_support_start("NEON")
        result = self.hasfunction(
            "int16x4_t acc = vdup_n_s16(0);",
            include="<arm_neon.h>",
        )
        self._print_support_end("NEON", result)
        return result


################################################################################
# Writing version control information to the module
################################################################################


def git_version():
    # Return the git revision as a string
    # copied from numpy setup.py
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ["SYSTEMROOT", "PATH"]:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env["LANGUAGE"] = "C"
        env["LANG"] = "C"
        env["LC_ALL"] = "C"
        out = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            env=env,
        ).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(["git", "rev-parse", "HEAD"])
        GIT_REVISION = out.strip().decode("ascii")
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def copy_same_file_pass(source, destination):
    try:
        # Try to copy
        shutil.copy(source, destination)
    except shutil.SameFileError:
        # Ignore the error
        pass


class StaticLibrary(Extension):
    def __init__(self, *args, **kwargs):
        self.export_include = kwargs.pop("export_include", [])
        Extension.__init__(self, *args, **kwargs)


class build_ext(_build_ext):
    def build_extension(self, ext):
        if isinstance(ext, StaticLibrary):
            self.build_static_extension(ext)
        else:
            _build_ext.build_extension(self, ext)

    def copy_extensions_to_source(self):
        _extensions = self.extensions
        self.extensions = [e for e in _extensions if not isinstance(e, StaticLibrary)]
        _build_ext.copy_extensions_to_source(self)
        self.extensions = _extensions

    def build_static_extension(self, ext):
        from distutils import log

        sources = ext.sources
        if sources is None or not isinstance(sources, (list, tuple)):
            raise DistutilsSetupError(
                f"in 'ext_modules' option (extension '{ext.name}'), "
                "'sources' must be present and must be "
                "a list of source filenames",
            )
        sources = list(sources)

        ext_path = self.get_ext_fullpath(ext.name)
        depends = sources + ext.depends
        if not (self.force or newer_group(depends, ext_path, "newer")):
            log.debug("skipping '%s' extension (up-to-date)", ext.name)
            return
        else:
            log.info("building '%s' extension", ext.name)

        extra_args = ext.extra_compile_args or []
        macros = ext.define_macros[:]
        for undef in ext.undef_macros:
            macros.append((undef,))
        objects = self.compiler.compile(
            sources,
            output_dir=self.build_temp,
            macros=macros,
            include_dirs=ext.include_dirs,
            debug=self.debug,
            extra_postargs=extra_args,
            depends=ext.depends,
        )
        self._built_objects = objects[:]
        if ext.extra_objects:
            objects.extend(ext.extra_objects)
        extra_args = ext.extra_link_args or []

        language = ext.language or self.compiler.detect_language(sources)

        libname = os.path.basename(ext_path).split(os.extsep)[0]
        output_dir = os.path.dirname(ext_path)

        if self.compiler.static_lib_format.startswith("lib") and libname.startswith(
            "lib",
        ):
            libname = libname[3:]

        # 1. copy to build directory
        # 1. copy to src tree for develop mode
        # 1. catch errors because we're running in editable mode
        import re

        src_tree_output_dir = re.match("build.*(mdtraj.*)", output_dir)

        try:
            # This might fail if no match in editable mode.
            src_tree_output_dir = src_tree_output_dir.group(1)
            if not os.path.exists(src_tree_output_dir):
                os.makedirs(src_tree_output_dir)
        except AttributeError:
            src_tree_output_dir = output_dir

        if not os.path.exists(output_dir):
            # necessary for windows
            os.makedirs(output_dir)

        assert os.path.isdir(src_tree_output_dir)

        self.compiler.create_static_lib(
            objects,
            output_libname=libname,
            output_dir=output_dir,
            target_lang=language,
        )

        lib_path = self.compiler.library_filename(libname, output_dir=output_dir)

        copy_same_file_pass(lib_path, src_tree_output_dir)

        for item in ext.export_include:
            copy_same_file_pass(item, src_tree_output_dir)
            copy_same_file_pass(item, output_dir)

    def get_ext_filename(self, ext_name):
        filename = _build_ext.get_ext_filename(self, ext_name)

        try:
            exts = [e for e in self.extensions if ext_name in {e.name, e.name.split(".")[-1]}]
            ext = exts[0]
            if isinstance(ext, StaticLibrary):
                if new_compiler().compiler_type == "msvc":
                    return filename.split(".")[0] + ".lib"
                else:
                    return filename.split(".")[0] + ".a"
        except Exception as e:  # noqa
            pass
        return filename


def parse_setuppy_commands():
    """Check the commands and respond appropriately.
    Return a boolean value for whether or not to run the build or not (avoid
    parsing Cython and template files if False).

    Adopted from scipy setup
    """
    args = sys.argv[1:]

    if not args:
        # User forgot to give an argument probably, let setuptools handle that.
        return True

    info_commands = [
        "--help-commands",
        "--name",
        "--version",
        "-V",
        "--fullname",
        "--author",
        "--author-email",
        "--maintainer",
        "--maintainer-email",
        "--contact",
        "--contact-email",
        "--url",
        "--license",
        "--description",
        "--long-description",
        "--platforms",
        "--classifiers",
        "--keywords",
        "--provides",
        "--requires",
        "--obsoletes",
    ]

    for command in info_commands:
        if command in args:
            return False

    # Note that 'alias', 'saveopts' and 'setopt' commands also seem to work
    # fine as they are, but are usually used together with one of the commands
    # below and not standalone.  Hence they're not added to good_commands.
    good_commands = (
        "develop",
        "sdist",
        "build",
        "build_ext",
        "build_py",
        "build_clib",
        "build_scripts",
        "bdist_wheel",
        "bdist_rpm",
        "bdist_wininst",
        "bdist_msi",
        "bdist_mpkg",
        "build_sphinx",
    )

    for command in good_commands:
        if command in args:
            return True

    # The following commands are supported, but we need to show more
    # useful messages to the user
    if "install" in args:
        return True

    if "--help" in args or "-h" in sys.argv[1]:
        return False

    # Commands that do more than print info, but also don't need Cython and
    # template parsing.
    other_commands = ["egg_info", "install_egg_info", "rotate"]
    for command in other_commands:
        if command in args:
            return False

    # If we got here, we didn't detect what setup.py command was given
    warnings.warn(
        "Unrecognized setuptools command ('{}'), proceeding with "
        "generating Cython sources and expanding templates".format(
            " ".join(sys.argv[1:]),
        ),
    )
    return True
