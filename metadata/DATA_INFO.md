## WHERE IS THE DATA

<!-- Data folder: `/nfs/team205/ed6/data/Fetal_immune/` -->

### Raw data matrices
- STARSOLO Mapping and cellbender outputs: `/lustre/scratch117/cellgen/team205/sharedData/ly5/cellbender/`
- Filtered `.h5ad` for each sample: `/nfs/team205/ed6/data/Fetal_immune/cellbender_raw/`
- Souporcell output: `...`

### Working anndata objects

Current data timestamp: `20210429`

- Pan Fetal atlas - raw counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_raw_count.{timestamp}.h5ad` 
- Pan Fetal atlas - log-normalized counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.h5ad` 
- Full obs table (collecting metadata and pre-existing annotations) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.full_obs.csv`
- var table (gene dispersion estimates for feature selection) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.var.csv`
- Pan Fetal atlas - log-normalized counts - with embeddings (`X_scVI` and `X_umap`) `/nfs/team205/ed6/data/Fetal_immune/scVI*.npy` 


