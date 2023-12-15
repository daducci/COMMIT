from commit.core import Evaluation, setup

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version
__version__ = version('dmri-commit')

__all__ = [Evaluation, setup, __version__]
