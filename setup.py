"""PEAK-QC a quality control tool for ATAC-seq data."""

from setuptools import setup
import re
import os


def find_version(f: str) -> str:
    """
    Get package version from version file.

    Parameters
    ----------
    f : str
        Path to version file.

    Returns
    -------
    str
        Version string.

    Raises
    ------
    RuntimeError
        If version string is missing.
    """
    version_file = open(f).read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="peakqc",
    description='Module for quality control of ATAC-seq data',
    version=find_version(os.path.join("peakqc", "_version.py")),
    license='MIT',
    packages=['peakqc'],
    python_requires='>=3.9',
    install_requires=["numpy",
                      "pandas",
                      "matplotlib",
                      "tqdm",
                      "beartype",
                      "matplotlib",
                      "scanpy>=1.9",
                      "pysam",
                      "scipy"])
