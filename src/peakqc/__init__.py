"""PEAKQC (peakqc)."""

# import with prefix _ to hide them
from ._version import __version__

# define what is exported in this module
__all__ = [
    "__version__"
]


def __dir__():
    """Return the defined submodules."""
    return __all__
