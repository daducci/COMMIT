import sys
from os import environ
import pyximport
from setuptools import Extension
from Cython.Build.Inline import _get_build_extension
from Cython.Build import cythonize
import os
from amico.models import StickZeppelinBall as _StickZeppelinBall, CylinderZeppelinBall as _CylinderZeppelinBall, VolumeFractions as _VolumeFractions

try:
    print('\n\nCOMMIT_WIP_MODEL\n\n')
    sys.path.append(environ["WIP_MODEL"])
    extension = Extension(name='commitwipmodels', language='c++', sources=[os.environ['WIP_MODEL'] + '/commitwipmodels.pyx'])
    build_extension = _get_build_extension()
    build_extension.extensions = cythonize([extension],
                                        include_path=[],
                                        quiet=False)
    build_extension.build_temp = os.environ['WIP_MODEL'] + '/build'
    build_extension.build_lib = os.environ['WIP_MODEL']

    build_extension.run()
    from commitwipmodels import *

except KeyError:
    print("KeyError: 'COMMIT_WIP_MODEL' not found in environ")
    pass
except ImportError:
    print("ImportError: 'commitwipmodels' not found in sys.modules")
    pass


class StickZeppelinBall(_StickZeppelinBall):
    """Simulate the response functions according to the Stick-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    pass


class CylinderZeppelinBall(_CylinderZeppelinBall):
    """Simulate the response functions according to the Cylinder-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    pass

class VolumeFractions(_VolumeFractions):
    """Simulate the response functions according to the VolumeFractions model.
    See the AMICO.model module for details.
    """
    pass
