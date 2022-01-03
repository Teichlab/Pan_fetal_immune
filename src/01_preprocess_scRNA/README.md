This folder contains scripts to merge scRNA-seq datasets in Pan Fetal Immune Atlas.

1. `preprocess_scRNA_uniform_old_labels.ipynb` - collapsing old annotations from publications or original authors 
2. `preprocess_scRNA_1_read_cellbender.py` - read `.h5` output from cellbender and run cell and gene filtering
3. `preprocess_scRNA_2_merge_cellbender.py` - merge filtered matrices for each sample in one anndata, remove doublets, save anndatas of raw counts and log-normalized counts matrices
4. `preprocess_scRNA_3_make_obs.py` - merge cell metadata and uniformed annotations and save table that becomes `adata.obs`
5. `preprocess_scRNA_4_make_var.py` - calculate highly variable genes and save table that becomes `adata.var`
6. `preprocess_scRNA_plot_metadata.ipynb` - calculate summary stats and make plots on metadata information
7. `souporcell_maternal_contams/` - folder containing scripts for detection of maternal contaminants through genotyping with [souporcell](https://github.com/wheaton5/souporcell) (generates [`Pan_fetal_immune/metadata/maternal_contaminants.csv`]())