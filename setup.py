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
    spline_smoothing = Extension(name='commit.spline_smoothing',
                     sources=['commit/spline_smoothing/spline_smoothing.pyx'],
                     include_dirs=['commit/trk2dictionary', 'commit/spline_smoothing/psimpl_v7_src'],
                     extra_compile_args=['-w'],
                     language='c++')
    bundle_o_graphy = Extension(name='commit.bundle_o_graphy',
                     sources=['commit/bundle_o_graphy.pyx'],
                     include_dirs=['commit/trk2dictionary', 'commit/spline_smoothing', 'commit/spline_smoothing/psimpl_v7_src'],
                     extra_compile_args=['-w'],
                     language='c++')

    return [trk2dictionary, core, proximals, spline_smoothing, bundle_o_graphy]

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
            self.parallel = 1
        print( f'Parallel compilation using {self.parallel} threads' )

        # Call original build_ext command
        build_ext.finalize_options(self)
        build_ext.run(self)

# import details from {package_name}/info.py
import sys
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
    install_requires=['setuptools>=46.1', 'Cython>=0.29', 'numpy>=1.12', 'scipy>=1.0', 'dipy>=1.0', 'dmri-amico>=1.3.2'],
    package_data={f'{package_name}.operator': ["*.*"]}
)
