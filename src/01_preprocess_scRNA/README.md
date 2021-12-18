_This folder contains scripts to merge datasets in Pan Fetal Immune Atlas._ 

1. `preprocess_scRNA_uniform_old_labels.ipynb` - collapsing old annotations from publications or original authors 
2. `preprocess_scRNA_1_read_cellbender.py` - reads `.h5` output from cellbender and runs cell and gene filtering
3. `preprocess_scRNA_2_merge_cellbender.py` - merges filtered matrices for each sample in one anndata, removes doublets, saves anndatas of raw counts and log-normalized counts matrices
4. `preprocess_scRNA_3_make_obs.py` - merges cell metadata and uniformed annotations and saves table that becomes `adata.obs`
5. `preprocess_scRNA_4_make_var.py` - calculates highly variable genes and saves table that becomes `adata.var`
6. `preprocess_scRNA_plot_metadata.ipynb` - calculate summary stats and make plots on metadata information