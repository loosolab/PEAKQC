# PEAKQC

Periodicity Evaluation As Key aspect of ATAC-seq Quality Control

A python tool for ATAC-seq quality control in single cells. 
On the bulk level quality control approaches rely on four key aspects: 

    - signal-to-noise ratio 
    - library complexity
    - mitochondrial DNA nuclear DNA ratio 
    - fragment length distribution 

Hereby relies PEAKQC on the evaluation of the fragment length distribution.
While on the bulk level the evaluation is done visually, it is not possible to do that on the single cell level.
PEAKQC solves this constraint with an convolution based algorithmic approach.


# Workflow

To execute the tool an anndata object and fragments, corresponding to the cells in the anndata have to be provided. The fragments can be either determined from a bamfile directly or by an fragments file in the bed format. If a fragments bedfile is available this is recommended to shorten the runtime.

![](/figures/PEAKQC_workflow.drawio.png)


# Installation

## PyPi
```
pip install peakqc
```
## From Source

### 1. Enviroment & Package Installation
1. Download the repository. This will download the repository to the current directory
```
git@gitlab.gwdg.de:loosolab/software/peakqc.git
```
2. Change the working directory to the newly created repository directory.
```
cd sc_framework
```
3. Install analysis environment. Note: using `mamba` is faster than `conda`, but this requires mamba to be installed.
```
mamba env create -f peakqc_env.yml
```
4. Activate the environment.
```
conda activate peakqc
```
5. Install PEAKQC into the enviroment.
```
pip install .
```

### 2. Package Installation
1. Download the repository. This will download the repository to the current directory
```
git@gitlab.gwdg.de:loosolab/software/peakqc.git
```
2. Change the working directory to the newly created repository directory.
```
cd sc_framework
```
3. Install PEAKQC into the enviroment.
```
pip install .
```

# Quickstart

Import the `add_fld_metrics()` function from the `fld_scoring` module.
Assure that the anndata.obs table contains the same barcodes as the provided fragment source. 
This includes that the barcodes are in the same format.

The anndata.obs column where the barcodes are stored can be specified with the `barcode_col` argument. When the barcodes are the index of the table leave the `barcode_col` argument to be `None`.

Next provide either a bedfile containing the fragments or a bamfile by the `fragments` argument. If available using a bedfile is recommended due to shorter runtime.

```
from peakqc.fld_scoring import *

adata = add_fld_metrics(adata=anndata,
                        fragments=fragments_bedfile,
                        barcode_col=None,
                        plot=True,
                        save_density=None,
                        save_overview=None,
                        sample=0,
                        n_threads=8)
```

