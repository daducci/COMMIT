from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy
import os

# os.environ["CC"] = "clang" 
# os.environ["CXX"] = "clang"

ext = Extension(
    name='spline_smoothing',
    sources=['spline_smoothing.pyx'],
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-w', '-std=c++11'],
    extra_link_args=[],
    language='c++',
)

setup(
    name='spline_smoothing',
    description='Apply spline filtering to each fiber in the input tractogram',
    version='1.0',
    cmdclass = {'build_ext':build_ext},
    ext_modules = [ ext ],
)
