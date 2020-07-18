from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import amico
import os
from os.path import join as pjoin

amico_version = amico.__version__.split('.')
amico_version = [int(version_val) for version_val in amico_version]
if amico_version[0] == 1 and amico_version[1] < 1:
    raise RuntimeError( 'COMMIT requires AMICO v1.1.0 or above. Current AMICO version is %s' % amico.__version__ )


# taken from npcuda
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

# Obtain the numpy include directory. This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# Try to locate CUDA
CUDA = locate_cuda()

if CUDA != None:
    # Run the customize_compiler
    class cuda_build_ext(build_ext):
        def build_extensions(self):
            customize_compiler_for_nvcc(self.compiler)
            build_ext.build_extensions(self)

    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    trk2dictionary_ext = Extension(
        name='commit.trk2dictionary',
        sources=['commit/trk2dictionary/trk2dictionary.pyx'],
        include_dirs=[numpy.get_include()],
        extra_compile_args= {
            'gcc': ['-w'],
            'nvcc': [
                '-arch=sm_30', '--ptxas-options=-v', '-c',
                '--compiler-options', "'-fPIC'"
                ]
            },
        extra_link_args=[],
        language='c++',
    )

    core_ext = Extension(
        name='commit.core',
        sources=['commit/core.pyx'],
        include_dirs=[numpy.get_include()],
        extra_compile_args= {
            'gcc': ['-w'],
            'nvcc': [
                '-arch=sm_30', '--ptxas-options=-v', '-c',
                '--compiler-options', "'-fPIC'"
                ]
            },
        extra_link_args=[],
        language='c++',
    )

    proximals_ext = Extension(
        name='commit.proximals',
        sources=['commit/proximals.pyx'],
        include_dirs=[numpy.get_include()],
        extra_compile_args= {
            'gcc': ['-w'],
            'nvcc': [
                '-arch=sm_30', '--ptxas-options=-v', '-c',
                '--compiler-options', "'-fPIC'"
                ]
            },
        extra_link_args=[],
        language='c++',
    )

    cudaoperator_ext = Extension(
        name='commit.cudaoperator',
        sources = ['commit/operator_withCUDA.cu', 'commit/cudaoperator.pyx'],
        library_dirs = [CUDA['lib64']],
        libraries = ['cudart'],
        language = 'c++',
        runtime_library_dirs = [CUDA['lib64']],
        # This syntax is specific to this build system
        # we're only going to use certain compiler args with nvcc
        # and not with gcc the implementation of this trick is in
        # customize_compiler()
        extra_compile_args= {
            'gcc': ['-w'],
            'nvcc': [
                '-arch=sm_30', '--ptxas-options=-v', '-c',
                '--compiler-options', "'-fPIC'"
                ]
            },
        include_dirs = [numpy_include, CUDA['include']]
    )

    setup(
        name='commit',
        version='1.4.0',
        description='Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)',
        author='Alessandro Daducci',
        author_email='alessandro.daducci@gmail.com',
        url='https://github.com/daducci/COMMIT',
        cmdclass = {'build_ext':cuda_build_ext},
        ext_modules = [ trk2dictionary_ext, core_ext, proximals_ext, cudaoperator_ext ],
        packages=['commit','commit.operator'],
        package_data={
            'commit.operator':["*.*"], # needed by pyximport to compile at runtime
        },
    )
else:
    print('Installing CPU version')

    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    ext1 = Extension(
        name='commit.trk2dictionary',
        sources=['commit/trk2dictionary/trk2dictionary.pyx'],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-w'],
        extra_link_args=[],
        language='c++',
    )

    ext2 = Extension(
        name='commit.core',
        sources=['commit/core.pyx'],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-w'],
        extra_link_args=[],
        language='c++',
    )

    ext3 = Extension(
        name='commit.proximals',
        sources=['commit/proximals.pyx'],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-w'],
        extra_link_args=[],
        language='c++',
    )

    setup(
        name='commit',
        version='1.3.0',
        description='Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)',
        author='Alessandro Daducci',
        author_email='alessandro.daducci@gmail.com',
        url='https://github.com/daducci/COMMIT',
        cmdclass = {'build_ext':build_ext},
        ext_modules = [ ext1, ext2, ext3 ],
        packages=['commit','commit.operator'],
        package_data={
            'commit.operator':["*.*"], # needed by pyximport to compile at runtime
        },
    )
