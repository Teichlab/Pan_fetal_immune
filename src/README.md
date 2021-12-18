This directory contains source code used in the preprocessing of data or analysis of the project. 

Each subfolders contains scripts and notebooks used for a step in the analysis. The numbering indicates the order in which analysis steps are executed. Files/folders that are not numbered and start with a underscore `_` are not part of the main analysis steps.  

## Contents
* `01_preprocess_scRNA`: scripts to collect raw data matrices in a single anndata object and collecting pre-existing annotations and metadata.
* `02_integration_scRNA`: scripts for data integration (with BBKNN or scVI), common embedding and clustering, splitting into lineage subsets
* `03_annotation`: scripts for cell type annotation based on marker gene expression
* `5_organ_signatures`: scripts for analysis of organ-specific cell type signatures (with LMM and factor analysis)
* `6_trajectory_inference`: scripts for trajectory inference on immune cells
* `7_differential_abundance`: scripts for differential abundance analysis in time with Milo
* `8_align_query`: scripts to map new data to panfetal references with scArches implementation in `scvi-tools`
* `_misc` miscellaneous analyses and scripts
* `utils` contains utility functions used in many analyses
