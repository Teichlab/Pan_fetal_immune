This directory contains scripts and notebooks used in the preprocessing and analysis of fetal immune dataset. 

## Contents

* `01_preprocess_scRNA`: scripts to collect raw data matrices in a single anndata object and collecting pre-existing annotations and metadata.
* `02_integration_scRNA`: scripts for data integration (with [scVI](https://scvi-tools.org/) or [BBKNN](https://github.com/Teichlab/bbknn)), common embedding and clustering, splitting into lineage subsets
* `03_annotation`: scripts for cell type annotation based on marker gene expression
* `04_milo_analysis`: scripts and notebooks for differential abundance analysis with [Milo](https://github.com/emdann/milopy) and differential expression analysis of variation across organs and gestation in lymphoid and myeloid compartments.
* `05_adult2fetal_mapping`: scripts for mapping of adult data from [Pan Immune Project](https://www.biorxiv.org/content/10.1101/2021.04.28.441762v2.full) to fetal reference with [scArches](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scarches_scvi_tools.html)
* `06_spatial`: scripts for spatial mapping of cell types in Visium spatial transcriptomics data with [cell2location](https://cell2location.readthedocs.io/en/latest/)
* `07_widespread_hematopoiesis`: scripts and notebooks for analysis of distribution of immune cell progenitors across organs, including analysis of spatial distribution of B cell progenitors and cell-cell interactions with [CellPhoneDB](https://github.com/Teichlab/cellphonedb).
* `08_B1`: scripts and notebooks for analysis of features of putative B1 cells in developing organs, including BCR analysis with [dandelion](https://github.com/zktuong/dandelion)
* `09_unconv_Tcells`: scripts and notebooks for characterization of unconventional T cells in developing organs, including TCR analysis with [scirpy](https://github.com/icbi-lab/scirpy)
* `_misc`: miscellaneous analyses and scripts
* `utils` contains utility functions and scripts used in multiple analyses
