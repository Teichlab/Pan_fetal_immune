## WHERE IS THE DATA

<!-- Data folder: `/nfs/team205/ed6/data/Fetal_immune/` -->

### Raw data matrices
- STARSOLO Mapping and cellbender outputs: `/lustre/scratch117/cellgen/team205/sharedData/ly5/cellbender/`
- Filtered `.h5ad` for each sample: `/nfs/team205/ed6/data/Fetal_immune/cellbender_raw/`
- Souporcell output: `...`

### Working anndata objects

Current data timestamp: `20210429`

#### Full atlas

- Pan Fetal atlas - raw counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_raw_count.{timestamp}.h5ad` 
- Pan Fetal atlas - log-normalized counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.h5ad` 
- Full obs table (collecting metadata and pre-existing annotations) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.full_obs.csv`
- var table (gene dispersion estimates for feature selection) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.var.csv`
- scVI latent embedding `/nfs/team205/ed6/data/Fetal_immune/scVI_outs/PAN.A01.v01.entire_data_raw_count.{t}.scVI_out.npy` 
- scVI model (to add new query datasets) `~/mount/gdrive/Pan_fetal/data4gpu_node/...`

#### Lineage subsets

Split IDs:
- `STROMA` (non-immune cells)
- `HSC_IMMUNE` (myeloid + lymphoid + erythroid + MK cells + progenitors)
- `MYELOID_LYMPHOID` (myeloid + lymphoid + progenitors)
- `MEM_PROGENITORS` (erythroid + MK cells + progenitors)
- `MYELOID` (myeloid + progenitors)
- `LYMPHOID` (lymphoid + progenitors)

Working anndata (log-norm expression + scVI embeddings + old annotations): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.{splitID}.embedding.h5ad` 
 
AnnData components:
- raw counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_raw_count.{timestamp}.{splitID}.h5ad` 
- var table (gene dispersion estimates for feature selection) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.{splitID}.var.csv`
- scVI latent embedding `/nfs/team205/ed6/data/Fetal_immune/scVI_outs/PAN.A01.v01.entire_data_raw_count.{t}.{s}.scVI_out.npy`  
- scVI model (to add new query datasets) `~/mount/gdrive/Pan_fetal/data4gpu_node/scvi_{splitID}_model`
 
