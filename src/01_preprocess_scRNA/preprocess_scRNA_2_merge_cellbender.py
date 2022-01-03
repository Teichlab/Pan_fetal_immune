### Merge samples in one anndata ###
## WARNING! This script requires a lot of RAM (> 75GB)
import os
import numpy as np
import scipy as scipy
import scanpy as sc
import pandas as pd
import pickle as pkl
from collections import defaultdict, Counter
import joblib as jl
import gc
from datetime import datetime

def timestamp():
    return datetime.now().strftime("20%y%m%d")

## Adapted from scjp module
def merge_matrix(ad, obskeys = None,use_raw = False,keep_only_mutual=False):
    '''merge matrix stored in ad
    ad: dictionary of anndata to merge
    obskeys: list to merge within anndata
    use_raw: if True, merge from .raw.X'''
    
    smp_list = list(ad.keys())
    obs_dict = defaultdict(list)
    obs_names = []
    
    for smp in smp_list:
        ad[smp].obs['name'] = smp
    
    if not obskeys:
        obskey_list = []
        obskeys = []
        for sample in smp_list:
            obskey_list.extend(list(ad[sample].obs.columns))
        for (obskey, number) in Counter(obskey_list).items():
            if number == len(smp_list):
                obskeys.append(obskey)
            else:
                if keep_only_mutual:
                    pass
                else:
                    for sample in smp_list:
                        if obskey not in ad[sample].obs.columns:
                            ad[sample].obs[obskey]='n/a'
                    obskeys.append(obskey)
                               
    for sample in smp_list:
        obs_names.extend(list(ad[sample].obs_names))
        for key in obskeys:   
            obs_dict[key].extend(list(ad[sample].obs[key]))
    
    from scipy.sparse import vstack
    if use_raw == True:
        stack = vstack([ad[x].raw.X for x in smp_list]) # stack data
        adata = sc.AnnData(stack, var = ad[smp_list[0]].raw.var)
    else:
        stack = vstack([ad[x].X for x in smp_list]) # stack data
        adata = sc.AnnData(stack, var = ad[smp_list[0]].var)
      
    
    adata.obs_names = obs_names
    print(len(adata))
    for obs_col in obs_dict:
        print(obs_col)
        adata.obs[obs_col] = obs_dict[obs_col]
    return adata

indir = '/nfs/team205/ed6/data/Fetal_immune/cellbender_raw/'
outdir = "/nfs/team205/ed6/data/Fetal_immune/"

## Read processed cellranger outputs
# (outputs of PFI_pp_1_read_cellbender.py)
print("Reading adatas...")
ad = {}
for filename in os.listdir(indir):
    ad[filename] = sc.read_h5ad(indir + filename)
    
## Merge
print("Merging adatas...")
merged_adata = merge_matrix(ad)

del ad
gc.collect()

## Remove doublets
print("Filtering doublets...")
merged_adata = merged_adata[merged_adata.obs["doublet_scores"] < 0.4]

## Write raw matrix
print("Saving raw matrix...")
merged_adata.write_h5ad(outdir + 'PAN.A01.v01.entire_data_raw_count.{t}.h5ad'.format(t=timestamp()))

## Normalize
print("Log-normalizing...")
sc.pp.normalize_per_cell(merged_adata, counts_per_cell_after=10e4)
sc.pp.log1p(merged_adata)

## Write normalized matrix
print("Saving normalized matrix...")
merged_adata.write_h5ad(outdir + 'PAN.A01.v01.entire_data_normalised_log.{t}.h5ad'.format(t=timestamp()))
