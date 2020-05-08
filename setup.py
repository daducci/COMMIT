from setuptools import setup
from setuptools.extension import Extension
import sys

try:
    from Cython.Build import cythonize
except ImportError:
    def cythonize(*args, **kwargs):
        from Cython.Build import cythonize
        return cythonize(*args, **kwargs)

try:
    from numpy import get_include
except ImportError:
    def get_include():
        from numpy import get_include
        return get_include()

with open('requirements.txt') as f:
    required_dependencies = f.read().splitlines()
    external_dependencies = []
    for dependency in required_dependencies:
        if dependency[0:2] == '-e':
            repo_name = dependency.split('=')[-1]
            repo_url = dependency[3:]
            external_dependencies.append('{} @ {}'.format(repo_name, repo_url))
        else:
            external_dependencies.append(dependency)

opts = dict(name='commit',
    version='1.3.7',
    description='Convex Optimization Modeling for Microstructure Informed Tractography (COMMIT)',
    author='Alessandro Daducci',
    author_email='alessandro.daducci@univr.it',
    url='https://github.com/daducci/COMMIT',
    packages=['commit','commit.operator'],
    setup_requires=['cython', 'numpy'],
    install_requires=external_dependencies,
    package_data={'commit.operator':["*.*"]},
)

def get_extensions():
    # Cython extension to create the sparse data structure from a tractogram
    # for the computation of matrix-vector multiplications
    ext1 = Extension(
        name='commit.trk2dictionary',
        sources=['commit/trk2dictionary/trk2dictionary.pyx'],
        include_dirs=[get_include()],
        extra_compile_args=['-w'],
        extra_link_args=[],
        language='c++',
    )

    ext2 = Extension(
        name='commit.core',
        sources=['commit/core.pyx'],
        include_dirs=[get_include()],
        extra_compile_args=['-w'],
        extra_link_args=[],
        language='c++',
    )

    ext3 = Extension(
        name='commit.proximals',
        sources=['commit/proximals.pyx'],
        include_dirs=[get_include()],
        extra_compile_args=['-w'],
        extra_link_args=[],
        language='c++',
    )
    return [ext1, ext2, ext3]


if __name__ == '__main__':
    setup(**opts)
    if sys.argv[1] != 'clean':
        opts['ext_modules'] = cythonize(get_extensions())
        setup(**opts)
