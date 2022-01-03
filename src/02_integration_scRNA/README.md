_This folder contains scripts for integration, common clustering and embedding of scRNA-seq datasets in Pan Fetal Immune Atlas._ 

1. `integration_scRNA_BBKNN.ipynb` - notebook for integration of full dataset with BBKNN (legacy version)
2. `integration_scRNA_scVI_fulldata.ipynb` - notebook for full dataset integration with [scVI](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html), exploration of results, qualitative comparison with BBKNN integration, splitting dataset into lineage subsets _(requires GPU)_
3. `integration_scRNA_run_scVI_subsets.py` - script to run scVI on dataset splits _(requires GPU)_
4. `integration_scRNA_cluster_scVI_subsets.py` - script for clustering and embedding on dataset splits after scVI integration
