from .core import Evaluation
__all__ = ['core','models','solvers','trk2dictionary']

from pkg_resources import get_distribution
__version__ = get_distribution('dmri-commit').version
