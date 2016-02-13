from __future__ import print_function, absolute_import
import os
import sys
import json
import shutil
import subprocess
import tempfile
from setuptools import Extension
from distutils.dep_util import newer_group
from distutils.errors import DistutilsExecError, DistutilsSetupError
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler, get_config_vars
from setuptools.command.build_ext import build_ext as _build_ext


################################################################################
# Detection of compiler capabilities
################################################################################

class CompilerDetection(object):
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

        self.msvc = cc.compiler_type == 'msvc'
        self._print_compiler_version(cc)

        if self.disable_openmp:
            self.openmp_enabled = False
        else:
            self.openmp_enabled, openmp_needs_gomp = self._detect_openmp()
        self.sse3_enabled = self._detect_sse3() if not self.msvc else True

        self.compiler_args_sse2 = ['-msse2'] if not self.msvc else ['/arch:SSE2']
        self.compiler_args_sse3 = ['-mssse3'] if (self.sse3_enabled and not self.msvc) else []
        self.compiler_args_warn = ['-Wno-unused-function', '-Wno-unreachable-code', '-Wno-sign-compare'] if not self.msvc else []

        if self.openmp_enabled:
            self.compiler_libraries_openmp = []

            if self.msvc:
                self.compiler_args_openmp = ['/openmp']
            else:
                self.compiler_args_openmp = ['-fopenmp']
                if openmp_needs_gomp:
                    self.compiler_libraries_openmp = ['gomp']
        else:
            self.compiler_libraries_openmp = []
            self.compiler_args_openmp = []

        if self.msvc:
            self.compiler_args_opt = ['/O2']
        else:
            self.compiler_args_opt = ['-O3', '-funroll-loops']
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
                cc.spawn([cc.compiler[0]] + ['-v'])
        except DistutilsExecError:
            pass

    def hasfunction(self, funcname, include=None, libraries=None, extra_postargs=None):
        # running in a separate subshell lets us prevent unwanted stdout/stderr
        part1 = '''
from __future__ import print_function
import os
import json
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler, get_config_vars

FUNCNAME = json.loads('%(funcname)s')
INCLUDE = json.loads('%(include)s')
LIBRARIES = json.loads('%(libraries)s')
EXTRA_POSTARGS = json.loads('%(extra_postargs)s')
        ''' % {
            'funcname': json.dumps(funcname),
            'include': json.dumps(include),
            'libraries': json.dumps(libraries or []),
            'extra_postargs': json.dumps(extra_postargs)}

        part2 = '''
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
        '''
        tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
        try:
            curdir = os.path.abspath(os.curdir)
            os.chdir(tmpdir)
            with open('script.py', 'w') as f:
                f.write(part1 + part2)
            proc = subprocess.Popen(
                [sys.executable, 'script.py'],
                stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            proc.communicate()
            status = proc.wait()
        finally:
            os.chdir(curdir)
            shutil.rmtree(tmpdir)

        return status == 0

    def _print_support_start(self, feature):
        print('Attempting to autodetect {0:6} support...'.format(feature), end=' ')

    def _print_support_end(self, feature, status):
        if status is True:
            print('Compiler supports {0}'.format(feature))
        else:
            print('Did not detect {0} support'.format(feature))

    def _detect_openmp(self):
        self._print_support_start('OpenMP')
        hasopenmp = self.hasfunction('omp_get_num_threads()', extra_postargs=['-fopenmp', '/openmp'])
        needs_gomp = hasopenmp
        if not hasopenmp:
            hasopenmp = self.hasfunction('omp_get_num_threads()', libraries=['gomp'])
            needs_gomp = hasopenmp
        self._print_support_end('OpenMP', hasopenmp)
        return hasopenmp, needs_gomp

    def _detect_sse3(self):
        "Does this compiler support SSE3 intrinsics?"
        self._print_support_start('SSE3')
        result = self.hasfunction('__m128 v; _mm_hadd_ps(v,v)',
                           include='<pmmintrin.h>',
                           extra_postargs=['-msse3'])
        self._print_support_end('SSE3', result)
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


def write_version_py(VERSION, ISRELEASED, filename='mdtraj/version.py'):
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


class StaticLibrary(Extension):
    def __init__(self, *args, **kwargs):
        self.export_include = kwargs.pop('export_include', [])
        Extension.__init__(self, *args, **kwargs)


class build_ext(_build_ext):
    def initialize_options(self):
        _build_ext.initialize_options(self)
        import pkg_resources
        dir = pkg_resources.resource_filename('numpy', 'core/include')
        self.include_dirs = [dir]

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
                  ("in 'ext_modules' option (extension '%s'), " +
                   "'sources' must be present and must be " +
                   "a list of source filenames") % ext.name)
        sources = list(sources)

        ext_path = self.get_ext_fullpath(ext.name)
        depends = sources + ext.depends
        if not (self.force or newer_group(depends, ext_path, 'newer')):
            log.debug("skipping '%s' extension (up-to-date)", ext.name)
            return
        else:
            log.info("building '%s' extension", ext.name)

        extra_args = ext.extra_compile_args or []
        macros = ext.define_macros[:]
        for undef in ext.undef_macros:
            macros.append((undef,))
        objects = self.compiler.compile(sources,
                                         output_dir=self.build_temp,
                                         macros=macros,
                                         include_dirs=ext.include_dirs,
                                         debug=self.debug,
                                         extra_postargs=extra_args,
                                         depends=ext.depends)
        self._built_objects = objects[:]
        if ext.extra_objects:
            objects.extend(ext.extra_objects)
        extra_args = ext.extra_link_args or []

        language = ext.language or self.compiler.detect_language(sources)

        libname = os.path.basename(ext_path).split(os.extsep)[0]
        output_dir = os.path.dirname(ext_path)

        if (self.compiler.static_lib_format.startswith('lib') and
            libname.startswith('lib')):
            libname = libname[3:]

        if not os.path.exists(output_dir):
            # necessary for windows
            os.makedirs(output_dir)

        self.compiler.create_static_lib(objects,
            output_libname=libname,
            output_dir=output_dir,
            target_lang=language)

        for item in ext.export_include:
            shutil.copy(item, output_dir)

    def get_ext_filename(self, ext_name):
        filename = _build_ext.get_ext_filename(self, ext_name)

        try:
            exts = [e for e in self.extensions if ext_name in {e.name, e.name.split('.')[-1]}]
            ext = exts[0]
            if isinstance(ext, StaticLibrary):
                if new_compiler().compiler_type == 'msvc':
                    return filename.split('.')[0] + '.lib'
                else:
                    return filename.split('.')[0] + '.a'
        except Exception as e:
            pass
        return filename
