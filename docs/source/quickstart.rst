Quickstart
==========

Example usage of `add_fld_metrics`:

.. code-block:: python

    from peakqc.fld_scoring import *

    adata = add_fld_metrics(adata=anndata,
                            fragments=fragments_bedfile,
                            barcode_col=None,
                            plot=True,
                            save_density=None,
                            save_overview=None,
                            sample=0,
                            n_threads=8)

Requirements:
- Barcodes must match between `anndata.obs` and fragment source
- Bed format is preferred for speed