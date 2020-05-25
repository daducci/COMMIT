from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


def get_extensions():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    ext1 = Extension(name='commit.trk2dictionary',
                     sources=['commit/trk2dictionary/trk2dictionary.pyx'],
                     extra_compile_args=['-w'],
                     extra_link_args=[],
                     language='c++')

    ext2 = Extension(name='commit.core',
                     sources=['commit/core.pyx'],
                     extra_compile_args=['-w'],
                     extra_link_args=[],
                     language='c++')

    ext3 = Extension(name='commit.proximals',
                     sources=['commit/proximals.pyx'],
                     extra_compile_args=['-w'],
                     extra_link_args=[],
                     language='c++')

    return [ext1, ext2, ext3]


class CustomBuildExtCommand(build_ext):
    """ build_ext command to use when numpy headers are needed. """

    def run(self):
        # Now that the requirements are installed, get everything from numpy
        from numpy import get_include
        from Cython.Build import cythonize

        # Add everything requires for build
        self.swig_opts = None
        self.include_dirs = [get_include()]
        self.extensions = cythonize(self.extensions)

        # Call original build_ext command
        build_ext.finalize_options(self)
        build_ext.run(self)


opts = dict(name='dmri-commit',
            version='1.3.8.1',
            description='Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)',
            long_description='Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)',
            author='Alessandro Daducci',
            author_email='alessandro.daducci@univr.it',
            url='https://github.com/daducci/COMMIT',
            packages=['commit', 'commit.operator'],
            cmdclass={'build_ext': CustomBuildExtCommand},
            ext_modules=get_extensions(),
            setup_requires=['Cython>=0.29.17', 'numpy>=1.18.4'],
            install_requires=['dmri-amico>=1.2.3', 'dipy>=1.1.0', 'Cython>=0.29.17', 'numpy>=1.18.4'],
            package_data={'commit.operator': ["*.*"]})

setup(**opts)

