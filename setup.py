from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
import sys

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
    """ build_ext command to use when numpy headers are needed. """
    def run(self):
        # Now that the requirements are installed, get everything from numpy
        from Cython.Build import cythonize
        from numpy import get_include
        from multiprocessing import cpu_count

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

# import details from {package_name}/info.py
sys.path.insert(0, f'./{package_name}/')
import info

# install the package
setup(
    name=info.NAME,
    version=info.VERSION,
    description=info.DESCRIPTION,
    long_description=info.LONG_DESCRIPTION,
    author=info.AUTHOR,
    author_email=info.AUTHOR_EMAIL,
    url=info.URL,
    license=info.LICENSE,
    packages=[f'{package_name}', f'{package_name}.operator'],
    cmdclass={'build_ext': CustomBuildExtCommand},
    ext_modules=get_extensions(),
    setup_requires=['Cython>=0.29', 'numpy>=1.12', 'wheel'],
    install_requires=['setuptools>=46.1', 'Cython>=0.29', 'numpy>=1.12', 'scipy>=1.0', 'dipy>=1.0', 'dmri-dicelib>=1.0.0', 'dmri-amico>=2.0.0'],
    package_data={f'{package_name}.operator': ["*.*"]}
)
