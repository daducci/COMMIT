import sys
from os import environ
from amico.models import StickZeppelinBall as _StickZeppelinBall, CylinderZeppelinBall as _CylinderZeppelinBall, VolumeFractions as _VolumeFractions

try:
    sys.path.append(environ['AMICO_WIP_MODELS'])
    from commitwipmodels import *
except KeyError:
    pass
except ImportError:
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
