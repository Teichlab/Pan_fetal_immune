## Prepare .h5ad files for cellxgene and submission
## Add metadata and annotations for the reference
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata
import scvelo as scv
import json

data_dir = '/nfs/team205/ed6/data/Fetal_immune/'
timestamp = '20210429'
output_dir = data_dir + 'cellxgene_h5ad_files/scRNA_data/'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

## Load full annotation table
anno_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.clean.csv".format(t=timestamp), index_col=0)

## Read barcodes of maternal cells
mat_barcodes = pd.read_csv("../../metadata/souporcell_results/maternal_barcodes.csv", index_col=0)
mat_barcodes = pd.Series([x.split("-1")[0] for x in mat_barcodes['x']])

## Load annotation groupings
with open('../../metadata/anno_groups.json', 'r') as json_file:
    anno_groups_dict = json.load(json_file)

all_splits = ['STROMA', 'HSC_IMMUNE', "MEM_PROGENITORS", "HSC_PROGENITORS","LYMPHOID", "MYELOID_V2","NKT"]

for split in all_splits:
    print("Processing {s}...".format(s=split))
    ## Read .h5ad file with embeddings
    adata = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.embedding.h5ad'.format(t=timestamp, s=split))
    adata.var_names_make_unique()

    ## Mark maternal contaminants
    adata.obs['is_maternal_contaminant'] = adata.obs_names.isin(mat_barcodes)

    ## Add fine cell type annotations
    adata.obs['anno_lvl_2_final_clean'] = np.nan
    adata.obs.loc[adata.obs_names.isin(anno_obs.index), 'anno_lvl_2_final_clean'] = anno_obs.loc[adata.obs_names[adata.obs_names.isin(anno_obs.index)]]['anno_lvl_2_final_clean']

    ## Rename some obs columns
    adata.obs['celltype_annotation'] = adata.obs['anno_lvl_2_final_clean'].copy()
    adata.obs['old_annotation_uniform'] = adata.obs['uniform_label_expanded_merged'].copy()

    keep_obs_cols = ['n_counts', 
                     'n_genes', 
                     'file', 
                     'mito', 
                     'doublet_scores',
                     'predicted_doublets', 
                     'old_annotation_uniform', 
                     'organ',
                     'Sort_id', 'age', 'method', 'donor', 'sex',
                     'Sample', 
                     'scvi_clusters', 
                     'is_maternal_contaminant',
                     'anno_lvl_2_final_clean', 
                     'celltype_annotation']

    adata.obs = adata.obs[keep_obs_cols].copy()
    adata.obs["donor"] = adata.obs["donor"].astype(str) ## Fix mixed type issue

    ## Clean uns
    adata.uns = {k:v for k,v in adata.uns.items() if not k.endswith("_colors")}

    ## Store the raw count matrix
    adata_raw = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.{s}.h5ad'.format(t=timestamp, s=split))
    adata.X = adata_raw.X.copy()

    ## Save info on vars used in scVI embedding 
    adata_var = pd.read_csv(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.var.csv'.format(t=timestamp, s=split), index_col=0)
    scvi_model_vars = pd.read_csv(output_dir + f'/scVI_models/scvi_{split}_model/var_names.csv', header=None)[0]
    adata_var['scvi_model_var'] = adata_var['GeneID'].isin(scvi_model_vars)
    adata.var = adata_var.copy()
    
    ## Save for data download
    adata.write_h5ad(output_dir + 'PAN.A01.v01.raw_count.{t}.{s}.embedding.h5ad'.format(t=timestamp, s=split), compression='gzip')
    
## Save objects 4 cellxgene
for split in all_splits:
    adata = sc.read_h5ad(output_dir + 'PAN.A01.v01.raw_count.{t}.{s}.embedding.h5ad'.format(t=timestamp, s=split))
    # log normalize
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10e4)
    sc.pp.log1p(adata)
    # Remove unnecessary components
    adata.uns = {}
    adata.obsp = {}
    del adata.obsm['X_scvi']
    adata.write_h5ad(output_dir + 'PAN.A01.v01.raw_count.{t}.{s}.embedding.cellxgene.h5ad'.format(t=timestamp, s=split))
