from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import amico

amico_version = amico.__version__.split('.')
amico_version = [int(version_val) for version_val in amico_version]
if amico_version[0] == 1 and amico_version[1] < 2:
    raise RuntimeError( 'COMMIT requires AMICO v1.2.0 or above. Current AMICO version is %s' % amico.__version__ )

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
    version='1.3.6',
    description='Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)',
    author='Alessandro Daducci',
    author_email='alessandro.daducci@univr.it',
    url='https://github.com/daducci/COMMIT',
    cmdclass = {'build_ext':build_ext},
    ext_modules = [ ext1, ext2, ext3 ],
    packages=['commit','commit.operator'],
    package_data={
        'commit.operator':["*.*"], # needed by pyximport to compile at runtime
    },
)
