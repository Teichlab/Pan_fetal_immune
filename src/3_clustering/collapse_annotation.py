### Merge and save manual annotations ### 

import os,sys
import numpy as np 
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata
import scipy
import matplotlib.pyplot as plt

## Load data
data_dir = '/nfs/team205/ed6/data/Fetal_immune/'
PFI_prefix = "PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118"

adata = sc.read_h5ad(data_dir + PFI_prefix + ".h5ad", backed="r")

### Load obs
new_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.wGut.full_obs.csv", index_col=0)
adata.obs = new_obs[new_obs.doublet_scores < 0.4]
cl_obs = pd.read_csv("/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118.clustering.obs.csv", index_col=0)
cl_obs.index = adata.obs_names
adata.obs = pd.concat([adata.obs, cl_obs.loc[adata.obs_names][["leiden_100", "leiden_150"]]],axis=1)

adata_obs = adata.obs
adata_var = adata.var

### Save info on splits in main anndata
# - Which cells belong to which split
# - which HVGS where used for each split
keep_splits = ["NKT", "B_CLEAN", "MYELOID_LYMPHOID", "STROMA", "MYELOID_CLEAN","LYMPHOID", "MEM_PROGENITORS"]

for s in keep_splits:
    split_file = "{prefix}.{split}.batchCorrected.h5ad".format(prefix = PFI_prefix, split = s)
    if split_file in os.listdir(data_dir):
        s_adata = sc.read_h5ad(data_dir + split_file, backed="r")
    else:
        print("No .h5ad found for split" + s)
    ## Store info on cell assignment to splits 
    adata_obs["isin_split_" + s] = adata_obs.index.isin(s_adata.obs_names)
    ## Store info on hvgs for this split
    adata_var["isHVG_in_split_" + s] = adata_var.index.isin(s_adata.var_names)
    del s_adata

adata_var.to_csv(data_dir + PFI_prefix  + ".full_var.csv")
    
### Load manual annotations
anno_dir = '/nfs/team205/ed6/bin/Pan_fetal_immune/manual_annotation/'

# keep_anno = [x.split(PFI_prefix)[1].split(".")[1] for x in os.listdir(anno_dir) if PFI_prefix in x]
keep_anno = ["LYMPHOID", "MYELOID", "MEM_PROGENITORS"]

for anno in keep_anno:
    adata_obs["anno_lvl_1_" + anno] = np.nan
    adata_obs["anno_lvl_2_" + anno] = np.nan
    anno_file = "{prefix}.{anno}.batchCorrected_annotation.csv".format(prefix = PFI_prefix, anno = anno)
    if anno_file in os.listdir(anno_dir):
        anno_df = pd.read_csv(anno_dir + anno_file, index_col=0)
    else:
        print("No .csv found for annotation" + anno)
    ## If there is no anno_lvl_2 replicate anno_lvl_1
    if not 'anno_lvl_1' in anno_df.columns:
        anno_df["anno_lvl_1"] = anno_df["anno_lvl_2"]
    if not 'anno_lvl_2' in anno_df.columns:
        anno_df["anno_lvl_2"] = anno_df["anno_lvl_1"]
    adata_obs.loc[anno_df.index,"anno_lvl_1_" + anno] = anno_df["anno_lvl_1"]
    adata_obs.loc[anno_df.index,"anno_lvl_2_"+ anno] = anno_df["anno_lvl_2"]

anno_cols = ["anno_lvl_1_" + a for a in keep_anno] + ["anno_lvl_2_" + a for a in keep_anno]

## Save info on missing annotations or duplicated annotations
adata_obs['is_annotated'] = adata_obs[anno_cols].notna().any(axis=1).astype("int")
adata_obs['is_uniquely_annotated'] = (adata_obs[anno_cols].notna().sum(axis=1)==2).astype("int")

## Save 
adata_obs.to_csv(data_dir + PFI_prefix  + ".full_obs.annotated.csv")
