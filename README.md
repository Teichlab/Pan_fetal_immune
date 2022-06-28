# Mapping the developing human immune system across organs

[![DOI](https://zenodo.org/badge/325343185.svg)](https://zenodo.org/badge/latestdoi/325343185)

Data processing and analysis scripts for fetal immune atlas (see our [paper](https://www.science.org/doi/10.1126/science.abo0510))

## Contents

* [`src`](https://github.com/emdann/Pan_fetal_immune/edit/master/src) contains scripts and notebooks used for data processing and analysis.
* [`metadata`](https://github.com/emdann/Pan_fetal_immune/edit/master/metadata): contains metadata relevant for sample and cell annotations, pointers to raw and processed data, color palettes and groupings.
* [`tutorials`](https://github.com/emdann/Pan_fetal_immune/edit/master/tutorials): contains tutorials for model re-use (details [here](https://github.com/Teichlab/Pan_fetal_immune/tree/metadata_curation#model-re-use-tutorials))

## Data and metadata information

Browse all processed datasets, models and annotations at [https://developmental.cellatlas.io/fetal-immune](https://developmental.cellatlas.io/fetal-immune).

### Processed gene expression data 

### Processed spatial transcriptomics data

### Processed scVDJ-seq data


### Raw sequencing data 
Raw sequencing libraries are deposited in ArrayExpress 
- scRNA-seq libraries: raw data for libraries generated for this study are deposited under accession [E-MTAB-11343](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11343/). Raw data for libraries published in previous studies can be found under accessions 
- []()      
- Visium libraries: [E-MTAB-11341](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11341/) 
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





