from importlib import metadata

from .core import Evaluation, setup

__all__ = ['core','models','solvers','trk2dictionary']

try:
    __version__ = metadata.version('dmri-commit')
except metadata.PackageNotFoundError:
    __version__ = 'not installed'