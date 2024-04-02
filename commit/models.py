import sys
from os import environ
import os
import pyximport
from setuptools import Extension
from Cython.Build.Inline import _get_build_extension
from Cython.Build import cythonize
import numpy as np
from amico.models import BaseModel, StickZeppelinBall as _StickZeppelinBall, CylinderZeppelinBall as _CylinderZeppelinBall
from amico.util import ERROR, set_verbose

try:
    sys.path.append(environ["WIP_MODEL"])
    extension = Extension(name='commitwipmodels',
                          language='c++',
                          sources=[os.environ['WIP_MODEL'] + '/commitwipmodels.pyx'],
                          extra_compile_args=['-w'])
    build_extension = _get_build_extension()
    build_extension.extensions = cythonize([extension],
                                        include_path=[],
                                        quiet=False)
    build_extension.build_temp = os.environ['WIP_MODEL'] + '/build'
    build_extension.build_lib = os.environ['WIP_MODEL']

    build_extension.run()
    from commitwipmodels import *

except ValueError:
    # check if .so exists
    path_files = os.listdir(environ["WIP_MODEL"])
    for f in path_files:
        if f.startswith('commitwipmodels') and f.endswith('.so'):
            from commitwipmodels import *
    else:
        pass

except KeyError:
    pass
except ImportError:
    pass


class StickZeppelinBall(_StickZeppelinBall):
    """Simulate the response functions according to the Stick-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        set_verbose(2)
        return super().resample(in_path, idx_out, Ylm_out, doMergeB0, ndirs)



class CylinderZeppelinBall(_CylinderZeppelinBall):
    """Simulate the response functions according to the Cylinder-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        set_verbose(2)
        return super().resample(in_path, idx_out, Ylm_out, doMergeB0, ndirs)

class VolumeFractions(BaseModel):
    """Implements a simple model where each compartment contributes only with
       its own volume fraction. This model has been created to test there
       ability to remove false positive fibers with COMMIT.
    """
    def __init__(self):
        self.id = 'VolumeFractions'
        self.name = 'Volume fractions'
        self.maps_name = []
        self.maps_descr = []

    def set(self):
        return

    def get_params(self):
        params = {}
        params['id'] = self.id
        params['name'] = self.name
        return params

    def set_solver(self):
        ERROR('Not implemented')

    def generate(self, out_path, aux, idx_in, idx_out, ndirs):
        return

    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        if doMergeB0:
            nS = 1 + self.scheme.dwi_count
            merge_idx = np.hstack((self.scheme.b0_idx[0], self.scheme.dwi_idx))
        else:
            nS = self.scheme.nS
            merge_idx = np.arange(nS)

        KERNELS = {}
        KERNELS['model'] = self.id
        KERNELS['wmr'] = np.ones((1, ndirs, nS), dtype=np.float32)
        KERNELS['wmh'] = np.ones((0, ndirs, nS), dtype=np.float32)
        KERNELS['iso'] = np.ones((0, nS), dtype=np.float32)
        return KERNELS

    def fit(self, evaluation):
        ERROR('Not implemented')