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
                     libraries=['stdc++'],
                     extra_compile_args=['-w', '-std=c++11'],
                     language='c++')
    core = Extension(name=f'{package_name}.core',
                     sources=[f'{package_name}/core.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')
    proximals = Extension(name=f'{package_name}.proximals',
                     sources=[f'{package_name}/proximals.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')
    operator = Extension(name=f'{package_name}.operator.operator',
                    sources=[f'{package_name}/operator/operator.pyx', f'{package_name}/operator/operator_c.c'],
                    extra_compile_args=['-w', '-O3', '-Ofast'],
                    language='c')
                 
    return [trk2dictionary, core, proximals, operator]

class CustomBuildExtCommand(build_ext):
    """ build_ext command to use when numpy headers are needed. """
    def run(self):
        # Now that the requirements are installed, get everything from numpy
        from multiprocessing import cpu_count

        from Cython.Build import cythonize
        from numpy import get_include

        # Add everything requires for build
        self.swig_opts = None
        self.include_dirs = [get_include()]
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
