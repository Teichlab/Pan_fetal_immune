This folder contains scripts to merge datasets in Pan Fetal Immune Atlas. See https://github.com/emdann/Pan_fetal_immune/blob/master/metadata/DATA_INFO.md for info on data location.

### Workflow description

#### Step 1: Read files `PFI_pp_1_read_cellbender.py`
Reads `.h5` output from cellbender and runs cell and gene filtering

#### Step 2: Merge in one anndata `PFI_pp_2_merge_cellbender.py`
Merges filtered matrices for each sample in one anndata, removes doublets, saves anndatas of raw counts and log-normalized counts matrices

#### Step 3: Make obs table `PFI_3_make_obs.py`
Merges metadata and uniformed annotations into one table that becomes the `adata.obs`
