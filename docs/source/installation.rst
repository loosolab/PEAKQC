Installation
============
PEAKQC can be either installed via the Python Package Index PyPI or from source.

Install from PyPi:
------------------

.. code-block:: bash

    pip install peakqc

Install from Source:
--------------------

.. code-block:: bash

    git clone git@gitlab.gwdg.de:loosolab/software/peakqc.git
    cd peakqc

Using conda/mamba:

.. code-block:: bash

    mamba env create -f peakqc_env.yml
    conda activate peakqc
    pip install .