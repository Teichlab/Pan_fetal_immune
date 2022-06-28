# Mapping the developing human immune system across organs

[![DOI](https://zenodo.org/badge/325343185.svg)](https://zenodo.org/badge/latestdoi/325343185)

Data processing and analysis scripts for fetal immune atlas (see our [paper](https://www.science.org/doi/10.1126/science.abo0510))

## Contents

* [`src`](https://github.com/emdann/Pan_fetal_immune/edit/master/src) contains scripts and notebooks used for data processing and analysis.
* [`metadata`](https://github.com/emdann/Pan_fetal_immune/edit/master/metadata): contains metadata relevant for sample and cell annotations, pointers to raw and processed data, color palettes and groupings.
* [`tutorials`](https://github.com/emdann/Pan_fetal_immune/edit/master/tutorials): contains tutorials for model re-use (details [here](https://github.com/Teichlab/Pan_fetal_immune/tree/metadata_curation#model-re-use-tutorials))

## Data and metadata information

Browse all processed datasets, models and annotations at [https://developmental.cellatlas.io/fetal-immune](https://developmental.cellatlas.io/fetal-immune).

### Primary tissue gene expression data  
For primary tissue samples, **gene expression count matrices** from scRNA-seq (10X Genomics 3' and 5') are available as [AnnData](https://anndata-tutorials.readthedocs.io/en/latest/getting-started.html) objects in `.h5ad` format for the [full dataset](https://developmental.cellatlas.io/fetal-immune#DatasetModalscRNA-seq1) (including stromal, immune and low quality cells), and 7 lineage subsets. In all objects `adata.X` stores gene counts corrected for background expression using [CellBender](https://github.com/broadinstitute/CellBender).

**Cell-level metadata**, including cell type annotations, are stored in `adata.obs` and can be downloaded separately [here](https://cellgeni.cog.sanger.ac.uk/developmentcellatlas/fetal-immune/PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv) (note: here cell type annotations are stored in column `anno_lvl_2_final_clean`). **Sample-level metadata**, including matching between scRNA-seq and scVDJ-seq libraries can be downloaded [here](). 

### _In vitro_ derived T cells gene expression data  
For in vitro derived T cells from Artificial Thymic Organoid protocols the **gene expression count matrix** from scRNA-seq (10X genomics 3') is available in `.h5ad` format [here](https://developmental.cellatlas.io/fetal-immune#DatasetModalscRNA-seq9). **Sample-level metadata** can be found [here](https://github.com/Teichlab/Pan_fetal_immune/blob/master/metadata/ATO_metadata_26102021.csv).

### Spatial transcriptomics data
Gene expression count matrices in `.h5ad` format and high-resolution images in `.tiff` format can be downloaded from the data portal (see section [Spatial Transcriptomics](https://developmental.cellatlas.io/fetal-immune)).  

#### Processed scVDJ-seq data


#### Raw sequencing data 
Raw sequencing libraries are deposited in ArrayExpress 
- scRNA-seq libraries: raw data for libraries generated for this study are deposited under accession [E-MTAB-11343](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11343/). Raw data for libraries published in previous studies can be found under accessions 
- []()      
- Visium libraries (incl. hi-res images): [E-MTAB-11341](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11341/) 
- scVDJ libraries: [E-MTAB-11388](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11388/)). 

  

## Model re-use tutorials

- [Mapping query data to fetal immune reference with scArches](https://nbviewer.org/github/Teichlab/Pan_fetal_immune/blob/master/tutorials/tutorial_query2reference_mapping.ipynb) 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/Pan_fetal_immune/blob/master/tutorials/tutorial_query2reference_mapping.ipynb)
- [Automatic cell type annotation from fetal immune reference with CellTypist](https://nbviewer.org/github/Teichlab/Pan_fetal_immune/blob/master/tutorials/tutorial_celltypist_fetal_immune.ipynb)[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/Pan_fetal_immune/blob/master/tutorials/tutorial_celltypist_fetal_immune.ipynb)

## Citation

If you use this data or code for your work, please cite

Suo C., Dann E., et al. (2022). Mapping the developing human immune system across organs. Science, 376(6597), https://doi.org/10.1126/science.abo0510

```
@article{suoMappingDevelopingHuman2022,
  title = {Mapping the Developing Human Immune System across Organs},
  author = {Suo, Chenqu and Dann, Emma and Goh, Issac and Jardine, Laura and Kleshchevnikov, Vitalii and Park, Jong-Eun and Botting, Rachel A. and Stephenson, Emily and Engelbert, Justin and Tuong, Zewen Kelvin and Polanski, Krzysztof and Yayon, Nadav and Xu, Chuan and Suchanek, Ondrej and Elmentaite, Rasa and Dom{\'i}nguez Conde, Cecilia and He, Peng and Pritchard, Sophie and Miah, Mohi and Moldovan, Corina and Steemers, Alexander S. and Mazin, Pavel and Prete, Martin and Horsfall, Dave and Marioni, John C. and Clatworthy, Menna R. and Haniffa, Muzlifah and Teichmann, Sarah A.},
  year = {2022},
  journal = {Science},
  volume = {376},
  number = {6597},
  publisher = {{American Association for the Advancement of Science}},
  doi = {10.1126/science.abo0510},
}
``` 

## Contact

For any questions, please post an [issue](https://github.com/emdann/Pan_fetal_immune/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc) in this repository or contact by email `ed6<at>sanger.ac.uk` or `cs42<at>sanger.ac.uk`.





