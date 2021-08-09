# WHERE IS THE DATA

## Raw data matrices
- STARSOLO Mapping and cellbender outputs: `/lustre/scratch117/cellgen/team205/sharedData/ly5/cellbender/`
- Filtered `.h5ad` for each sample: `/nfs/team205/ed6/data/Fetal_immune/cellbender_raw/`
- Souporcell output: `...`

### Working anndata objects

Current data timestamp: `20210429`

### Full atlas

- Pan Fetal atlas - raw counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_raw_count.{timestamp}.h5ad` 
- Pan Fetal atlas - log-normalized counts (doublets excluded): `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.h5ad` 
- Full obs table (collecting metadata and pre-existing annotations) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.full_obs.csv`
- var table (gene dispersion estimates for feature selection) `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.var.csv`
- scVI latent embedding `/nfs/team205/ed6/data/Fetal_immune/scVI_outs/PAN.A01.v01.entire_data_raw_count.{t}.scVI_out.npy` 
- scVI model (to add new query datasets) `~/mount/gdrive/Pan_fetal/data4gpu_node/...`
- Full obs table (collecting metadata and pre-existing annotations) **with manual annotations** `/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.{timestamp}.full_obs.annotated.clean.csv`

### Lineage subsets

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
- scVI latent embedding `/nfs/team205/ed6/data/Fetal_immune/scVI_outs/PAN.A01.v01.entire_data_raw_count.{timestamp}.{splitID}.scVI_out.npy`  
- scVI model (to add new query datasets) `~/mount/gdrive/Pan_fetal/data4gpu_node/scvi_{splitID}_model`

## Spatial data

- raw data: `/lustre/scratch117/cellgen/team205/sharedData/ed6/visium-noimage-hack` (mapped without image with a hack by `ktpolanski`)
- pre-processed anndata objects: `/nfs/team205/ed6/data/Fetal_immune/Visium/*.h5ad` 

## Analysis outputs

### Mapping of PIP data to Pan-fetal atlas

- AnnData of adult cells subsets `/nfs/team205/ed6/data/Fetal_immune/panimmune_{adult_subset}_query.h5ad`
- AnnData of adult cells mapped to fetal cells (with common embedding + score for similarity between query and reference cells in `adata.obs["mnn_sim_ratio"]`)`/nfs/team205/ed6/data/Fetal_immune/panimmune_{adult_subset}_query.mapped2{panfetal_subset}.withReference.h5ad` (where `adult_subset` denotes the subset of cells taken from adult, `panfetal_subset` the fetal subset)

### Differential abundance analysis with Milo

- dataframe of log-Fold Changes and FDR corrected P-vals by organ (as used in the beeswarm plots) `/nfs/team205/ed6/data/Fetal_immune/milo_outs/{split_ID}/milo_beeswarm_plot_data.{split_ID}.csv`
- AnnData objects with Milo specific components (as generated with `emdann/milopy` package) `/nfs/team205/ed6/data/Fetal_immune/milo_outs/{split_ID}/adata4milo.{split_ID}.csv`. This includes:
    - `adata.obsm["nhoods"]` a binary matrix matching cells (rows) to neighbourhoods (columns)
    - `adata.uns['nhood_adata']` an AnnData object that is used for differential abundance testing. Hhere `obs` are neighbourhoods, `vars` are samples and `adata.X` are counts of cells from each sample in each neighbourhood. This object is also separately saved in `/nfs/team205/ed6/data/Fetal_immune/milo_outs/{split_ID}/milo_nhood_adata.{split_ID}.h5ad`
