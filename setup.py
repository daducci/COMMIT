import os
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# name of the package
package_name = 'commit'

def get_extensions():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    trk2dictionary = Extension(name=f'{package_name}.trk2dictionary',
                     sources=[f'{package_name}/trk2dictionary/trk2dictionary.pyx'],
                     libraries=[] if sys.platform == 'win32' else ['stdc++'],
                     extra_compile_args=[] if sys.platform == 'win32' else ['-w', '-std=c++11'],
                     language='c++')
    core = Extension(name=f'{package_name}.core',
                     sources=[f'{package_name}/core.pyx'],
                     extra_compile_args=[] if sys.platform == 'win32' else ['-w', '-std=c++11'],
                     language='c++')
    proximals = Extension(name=f'{package_name}.proximals',
                     sources=[f'{package_name}/proximals.pyx'],
                     extra_compile_args=[] if sys.platform == 'win32' else ['-w', '-std=c++11'],
                     language='c++')
    models = Extension(name=f'{package_name}.models',
                     sources=[f'{package_name}/models.pyx'],
                     extra_compile_args=[] if sys.platform == 'win32' else ['-w', '-std=c++11'],
                     language='c++')    
    # NOTE: Windows requires the pthread-win32 static library to compile the operator extension
    #       The library can be downloaded from https://github.com/GerHobbelt/pthread-win32
    #       The PTHREAD_WIN_INCLUDE and PTHREAD_WIN_LIB environment variables must be set to the include and lib directories
    if sys.platform == 'win32':
        try:
            pthread_win_include = os.environ['PTHREAD_WIN_INCLUDE']
            pthread_win_lib = os.environ['PTHREAD_WIN_LIB']
        except KeyError:
            raise RuntimeError('PTHREAD_WIN_INCLUDE and PTHREAD_WIN_LIB must be set')
    operator = Extension(name=f'{package_name}.operator.operator',
                    sources=[f'{package_name}/operator/operator.pyx', f'{package_name}/operator/operator_c.c'],
                    include_dirs=[pthread_win_include] if sys.platform == 'win32' else [],
                    libraries=['pthread'] if sys.platform == 'win32' else [],
                    library_dirs=[pthread_win_lib] if sys.platform == 'win32' else [],
                    extra_compile_args=['/fp:fast', '/DHAVE_STRUCT_TIMESPEC'] if sys.platform == 'win32' else ['-w', '-O3', '-Ofast'],
                    language='c')
                 
    return [trk2dictionary, core, proximals, operator, models]

class CustomBuildExtCommand(build_ext):
    """ build_ext command to use when numpy headers are needed. """
    def run(self):
        # Now that the requirements are installed, get everything from numpy
        from multiprocessing import cpu_count

        from Cython.Build import cythonize
        from numpy import get_include

        # Add everything requires for build
        self.swig_opts = None
        self.include_dirs.extend([get_include()])
        self.distribution.ext_modules[:] = cythonize( self.distribution.ext_modules, build_dir='build' )

        # if not specified via '-j N' option, set compilation using max number of cores
        if self.parallel is None:
            self.parallel = cpu_count()
        print( f'Parallel compilation using {self.parallel} threads' )

        # Call original build_ext command
        build_ext.finalize_options(self)
        build_ext.run(self)

# generate the operator_c.c file
sys.path.insert(0, os.path.dirname(__file__))
from setup_operator import write_operator_c_file
write_operator_c_file()

# create the 'build' directory
if not os.path.exists('build'):
    os.makedirs('build')

# install the package
setup(
    cmdclass={'build_ext': CustomBuildExtCommand},
    ext_modules=get_extensions()
)
