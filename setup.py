""""PEAK-QC a quality control tool for ATAC-seq data"""

from setuptools import setup
from setuptools import find_namespace_packages
import re
import os
import glob

# Find all packages in sctoolbox
packages = find_namespace_packages("peak-qc")
packages = ["peak-qc." + package for package in packages]

# find top level scripts
#modules = glob.glob("peak-qc/*.py")
#modules = [m.replace("/", ".")[:-3] for m in modules if not m.endswith("__init__.py")]

def find_version(f: str) -> str:
    """
    Get package version from file.

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
    name="peak-qc",
    description='Module for quality control of ATAC-seq data',
    version=find_version(os.path.join("peak-qc", "_version.py")),
    license='MIT',
    packages=packages,
    python_requires='>=3.9',
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        ]
)