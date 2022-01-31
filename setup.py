from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

def get_extensions():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    trk2dictionary = Extension(name='commit.trk2dictionary',
                     sources=['commit/trk2dictionary/trk2dictionary.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')
    core = Extension(name='commit.core',
                     sources=['commit/core.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')
    proximals = Extension(name='commit.proximals',
                     sources=['commit/proximals.pyx'],
                     extra_compile_args=['-w'],
                     language='c++')
    return [trk2dictionary, core, proximals]

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

# import details from commit/info.py
import sys
sys.path.insert(0, './commit/')
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
    packages=['commit', 'commit.operator'],
    cmdclass={'build_ext': CustomBuildExtCommand},
    ext_modules=get_extensions(),
    setup_requires=['Cython>=0.29', 'numpy>=1.12'],
    install_requires=['wheel', 'setuptools>=46.1', 'Cython>=0.29', 'numpy>=1.12', 'scipy>=1.0', 'dipy>=1.0', 'dmri-amico>=1.3.2'],
    package_data={'commit.operator': ["*.*"]}
)
