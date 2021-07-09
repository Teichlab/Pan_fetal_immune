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
timestamp = '20210429'

PFI_prefix = "PAN.A01.v01.entire_data_normalised_log"

### Load obs
adata_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.csv".format(t=timestamp), index_col=0)

# ### Save info on splits in main anndata
# # - Which cells belong to which split
# # - which HVGS where used for each split
# keep_splits = ["NKT", "B_CLEAN", "MYELOID_LYMPHOID", "STROMA", "MYELOID_CLEAN","LYMPHOID", "MEM_PROGENITORS"]
    
### Load manual annotations
anno_dir = '/nfs/team205/ed6/bin/Pan_fetal_immune/metadata/manual_annotation/'

keep_anno = ["LYMPHOID", 
             "MYELOID", 
             "MEM_PROGENITORS"
            ]

for anno in keep_anno:
    adata_obs["anno_lvl_1_" + anno] = np.nan
    adata_obs["anno_lvl_2_" + anno] = np.nan
    anno_file = 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.csv'.format(t=timestamp, s=anno)
    if anno_file in os.listdir(anno_dir):
        anno_df = pd.read_csv(anno_dir + anno_file, index_col=0)
            ## If there is no anno_lvl_2 replicate anno_lvl_1
        if not 'anno_lvl_1' in anno_df.columns:
            anno_df["anno_lvl_1"] = anno_df["anno_lvl_2"]
        if not 'anno_lvl_2' in anno_df.columns:
            anno_df["anno_lvl_2"] = anno_df["anno_lvl_1"]
        adata_obs.loc[anno_df.index,"anno_lvl_1_" + anno] = anno_df["anno_lvl_1"]
        adata_obs.loc[anno_df.index,"anno_lvl_2_"+ anno] = anno_df["anno_lvl_2"]
    else:
        print("No .csv found for annotation " + anno)
        
anno_cols = ["anno_lvl_1_" + a for a in keep_anno] + ["anno_lvl_2_" + a for a in keep_anno]

## Save info on missing annotations or duplicated annotations
adata_obs['is_annotated'] = adata_obs[anno_cols].notna().any(axis=1).astype("int")
adata_obs['is_uniquely_annotated'] = (adata_obs[anno_cols].notna().sum(axis=1)==2).astype("int")

## Save 
adata_obs.to_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.csv".format(t=timestamp), index_col=0)