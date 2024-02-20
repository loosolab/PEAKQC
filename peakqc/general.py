"""General utility functions for the peakqc package."""

import importlib
from beartype import beartype
from beartype.typing import Any, Literal
import pysam


@beartype
def _is_gz_file(filepath: str) -> bool:
    """
    Check wheather file is a compressed .gz file.

    Parameters
    ----------
    filepath : str
        Path to file.

    Returns
    -------
    bool
        True if the file is a compressed .gz file.
    """

    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


@beartype
def check_module(module: str) -> None:
    """
    Check if <module> can be imported without error.

    Parameters
    ----------
    module : str
        Name of the module to check.

    Raises
    ------
    ImportError
        If the module is not available for import.
    """

    error = 0
    try:
        importlib.import_module(module)
    except ModuleNotFoundError:
        error = 1
    except Exception:
        raise  # unexpected error loading module

    # Write out error if module was not found
    if error == 1:
        s = f"ERROR: Could not find the '{module}' module on path, but the module is needed for this functionality. Please install this package to proceed."
        raise ImportError(s)


@beartype
def open_bam(file: str,
             mode: str,
             verbosity: Literal[0, 1, 2, 3] = 3, **kwargs: Any) -> pysam.AlignmentFile:
    """
    Open bam file with pysam.AlignmentFile. On a specific verbosity level.

    Parameters
    ----------
    file : str
        Path to bam file.
    mode : str
        Mode to open the file in. See pysam.AlignmentFile
    verbosity : Literal[0, 1, 2, 3], default 3
        Set verbosity level. Verbosity level 0 for no messages.
    **kwargs : Any
        Forwarded to pysam.AlignmentFile

    Returns
    -------
    pysam.AlignmentFile
        Object to work on SAM/BAM files.
    """

    # save verbosity, then set temporary one
    former_verbosity = pysam.get_verbosity()
    pysam.set_verbosity(verbosity)

    # open file
    handle = pysam.AlignmentFile(file, mode, **kwargs)

    # return to former verbosity
    pysam.set_verbosity(former_verbosity)

    return handle
