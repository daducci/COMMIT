from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


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


class CustomBuildExtCommand(build_ext):
    """ build_ext command to use when numpy headers are needed. """

    def run(self):
        # Now that the requirements are installed, get everything from numpy
        from Cython.Build import cythonize
        from numpy import get_include

        # Add everything requires for build
        self.swig_opts = None
        self.include_dirs = [get_include()]
        self.distribution.ext_modules[:] = cythonize(
                    self.distribution.ext_modules)

        # Call original build_ext command
        build_ext.finalize_options(self)
        build_ext.run(self)


description = 'Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)'
opts = dict(name='dmri-commit',
            version='1.5.1',
            description=description,
            long_description=description,
            author='Alessandro Daducci',
            author_email='alessandro.daducci@univr.it',
            url='https://github.com/daducci/COMMIT',
            license='BSD license',
            packages=['commit', 'commit.operator'],
            cmdclass={'build_ext': CustomBuildExtCommand},
            ext_modules=get_extensions(),
            setup_requires=['Cython>=0.29', 'numpy>=1.12'],
            install_requires=['Cython>=0.29', 'dmri-amico>=1.2.6', 'dipy>=1.0', 'numpy>=1.12'],
            package_data={'commit.operator': ["*.*"]})

setup(**opts)
