from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
import sys
import os

# name of the package
package_name = 'commit'

def get_extensions():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    libraries = []
    extra_compile_args = []
    language = 'c++'
    if sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        libraries = ['stdc++']
        extra_compile_args.extend(['-w', '-std=c++14'])
    if sys.platform.startswith('win32'):
        extra_compile_args.append('/std:c++14')
    
    trk2dictionary = Extension(name=f'{package_name}.trk2dictionary',
                    sources=[f'{package_name}/trk2dictionary/trk2dictionary.pyx'],
                    libraries=libraries,
                    extra_compile_args=extra_compile_args,
                    language=language)
    core = Extension(name=f'{package_name}.core',
                    sources=[f'{package_name}/core.pyx'],
                    extra_compile_args=extra_compile_args,
                    language=language)
    proximals = Extension(name=f'{package_name}.proximals',
                    sources=[f'{package_name}/proximals.pyx'],
                    extra_compile_args=extra_compile_args,
                    language=language)
        
    return [trk2dictionary, core, proximals]

class CustomBuildExtCommand(build_ext):
    def run(self):
        # Now that the requirements are installed, get everything from numpy
        from Cython.Build import cythonize
        from os import cpu_count

        # Add everything requires for build
        self.swig_opts = None
        self.distribution.ext_modules[:] = cythonize( self.distribution.ext_modules, build_dir='build' )

        # if not specified via '-j N' option, set compilation using max number of cores
        if self.parallel is None:
            self.parallel = cpu_count()
        print( f'Parallel compilation using {self.parallel} threads' )

        # Call original build_ext command
        build_ext.finalize_options(self)
        build_ext.run(self)

# create the 'build' directory
if not os.path.exists('build'):
    os.makedirs('build')

# install the package
setup(
    cmdclass={'build_ext': CustomBuildExtCommand},
    ext_modules=get_extensions()
)
