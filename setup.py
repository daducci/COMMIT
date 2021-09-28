from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
import os
from os.path import join as pjoin

# taken from https://github.com/rmcgibbo/npcuda-example/blob/master/cython/setup.py
def find_in_path(name, path):
    """Find a file in a search path"""

    # Adapted fom http://code.activestate.com/recipes/52224
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None

def locate_cuda():
    """Locate the CUDA environment on the system
    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.
    Starts by looking for the CUDAHOME env variable. If not found,
    everything is based on finding 'nvcc' in the PATH.
    """

    # First check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # Otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            return None
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home': home, 'nvcc': nvcc,
                  'include': pjoin(home, 'include'),
                  'lib64': pjoin(home, 'lib64')}
    for k, v in iter(cudaconfig.items()):
        if not os.path.exists(v):
            return None

    return cudaconfig

def customize_compiler_for_nvcc(self):
    """Inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.
    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on.
    """

    # Tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # Save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # Now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1
            # translated from the extra_compile_args in the Extension class
            print(type(extra_postargs))
            print(extra_postargs)
            postargs = extra_postargs['nvcc']
        else:
            print(type(extra_postargs))
            print(extra_postargs)
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # Reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # Inject our redefined _compile method into the class
    self._compile = _compile

# Locate CUDA
CUDA = locate_cuda()

def get_extensions():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    ext1 = Extension(name='commit.trk2dictionary',
                     sources=['commit/trk2dictionary/trk2dictionary.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')

    ext2 = Extension(name='commit.core',
                     sources=['commit/core.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')

    ext3 = Extension(name='commit.proximals',
                     sources=['commit/proximals.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')

    return [ext1, ext2, ext3]

def get_extensions_with_cuda():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications

    ext1 = Extension(name='commit.trk2dictionary',
                     sources=['commit/trk2dictionary/trk2dictionary.pyx'],
                     extra_compile_args= {'gcc':  ['-w'],
                                          'nvcc': ['-arch=sm_50', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                     extra_link_args=[],
                     language='c++')

    ext2 = Extension(name='commit.core',
                     sources=['commit/core.pyx'],
                     extra_compile_args= {'gcc':  ['-w'],
                                          'nvcc': ['-arch=sm_50', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                     extra_link_args=[],
                     language='c++')

    ext3 = Extension(name='commit.proximals',
                      sources=['commit/proximals.pyx'],
                      extra_compile_args= {'gcc':  ['-w'],
                                           'nvcc': ['-arch=sm_50', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                      extra_link_args=[],
                      language='c++')

    ext4 = Extension(name='commit.cudaoperator.operator',
                     sources = ['commit/cudaoperator/operator_withCUDA.cu', 'commit/cudaoperator/operator.pyx'],
                     extra_compile_args= {'gcc':  ['-w'],
                                          'nvcc': ['-arch=sm_50', '--ptxas-options=-v', '-c', '--compiler-options', "'-fPIC'"]},
                     language = 'c++',
                     library_dirs = [CUDA['lib64']],
                     libraries = ['cudart'],
                     runtime_library_dirs = [CUDA['lib64']])

    return [ext1, ext2, ext3, ext4]

if CUDA == None:
    extensions = get_extensions()
else:
    extensions = get_extensions_with_cuda()

if CUDA == None:
    class CustomBuildExtCommand(build_ext):
        """ build_ext command to use when numpy headers are needed. """

        def run(self):
            # Now that the requirements are installed, get everything from numpy
            from Cython.Build import cythonize
            from numpy import get_include
            
            # Add everything requires for build
            self.swig_opts = None
            self.include_dirs = [get_include()]
            self.distribution.ext_modules[:] = cythonize(self.distribution.ext_modules)

            # Call original build_ext command
            build_ext.finalize_options(self)
            build_ext.run(self)
else:
    class CustomBuildExtCommand(build_ext):
        """ build_ext command to use when numpy headers are needed. """

        def build_extensions(self):
            customize_compiler_for_nvcc(self.compiler)
            build_ext.build_extensions(self)

        def run(self):
            # Now that the requirements are installed, get everything from numpy
            from Cython.Build import cythonize
            from numpy import get_include
            
            # Add everything requires for build
            self.swig_opts = None
            self.include_dirs = [get_include(), CUDA['include'], 'commit/cudaoperator']
            self.distribution.ext_modules[:] = cythonize(self.distribution.ext_modules)

            # Call original build_ext command
            build_ext.finalize_options(self)
            build_ext.run(self)

description = 'Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)'
opts = dict(name='dmri-commit',
            version='1.5.0',
            description=description,
            long_description=description,
            author='Alessandro Daducci',
            author_email='alessandro.daducci@univr.it',
            url='https://github.com/daducci/COMMIT',
            license='BSD license',
            packages=['commit', 'commit.operator'],
            cmdclass={'build_ext': CustomBuildExtCommand},
            ext_modules=extensions,
            setup_requires=['Cython>=0.29', 'numpy>=1.12'],
            install_requires=['Cython>=0.29', 'dmri-amico>=1.2.6', 'dipy>=1.0', 'numpy>=1.12'],
            package_data={'commit.operator': ["*.*"]})

setup(**opts)