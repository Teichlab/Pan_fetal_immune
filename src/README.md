This directory contains source code used in the preprocessing of data or analysis of the project. 

Each subfolders contains scripts and notebooks used for a step in the analysis. The numbering indicates the order in which analysis steps are executed. Files/folders that are not numbered and start with a underscore `_` are not part of the main analysis steps.  

## Contents
* `01_preprocess_scRNA`: scripts to collect raw data matrices in a single anndata object and collecting pre-existing annotations and metadata.
* `02_integration_scRNA`: scripts for data integration (with BBKNN or scVI), common embedding and clustering, splitting into lineage subsets
* `03_annotation`: scripts for cell type annotation based on marker gene expression
* `04_milo_analysis`: scripts and notebooks for differential abundance and differential expression analysis of variation across organs and gestation in lymphoid and myeloid compartments 
* `05_adult2fetal_mapping`: scripts for mapping of adult data from [Pan Immune Project](https://www.biorxiv.org/content/10.1101/2021.04.28.441762v2.full) to fetal reference with [scArches](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scarches_scvi_tools.html)


* `_misc` miscellaneous analyses and scripts
* `utils` contains utility functions used in many analyses
