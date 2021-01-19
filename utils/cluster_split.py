### Clustering splits 4 annotation
import os,sys
import numpy as np 
import pandas as pd
import scanpy as sc
from datetime import datetime
import time

## import utils
cwd = '../utils/'
sys.path.append(cwd)

import genes
import panfetal_utils

import argparse
parser = argparse.ArgumentParser()
# parser.add_argument("split_adata", help="AnnData object for PFI atlas subset")
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma...)")
args = parser.parse_args()

## Load split data
save_path = '/nfs/team205/ed6/data/Fetal_immune/'
suffix = 'PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118'
s = args.split_name

adata_name = save_path + "{}.{}.h5ad".format(suffix, s)
sdata = sc.read_h5ad(adata_name)

sdata.obs['bbk'] = [x+y+z for x,y,z in zip(sdata.obs['donor'],sdata.obs['method'],sdata.obs['organ'])]

## Add updated obs
new_obs = pd.read_csv("/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.full_obs.csv", index_col=0)
sdata.obs = new_obs.loc[sdata.obs_names]

## Preprocess (feature selection, scale)
sdata_pp = panfetal_utils.pfi_preprocess(sdata, how="")

## Ridge regression
start=time.time()
panfetal_utils.ridge_regression(sdata_pp, batch_key=['method','donor','organ'],confounder_key=['uniform_label'])
v2_time = time.time() - start
print("Ridge regression runtime: ", str(v2_time))

sdata_pp.X = sdata_pp.layers['X_remain']
del sdata_pp.layers['X_remain']
del sdata_pp.layers['X_explained']

## BBKNN and clustering 
# Remove samples with less than 20 cells in split
sdata_pp_bbknn = sdata_pp[sdata_pp.obs['bbk'].isin(sdata_pp.obs["bbk"].value_counts().index[sdata_pp.obs["bbk"].value_counts() > 20])]
panfetal_utils.pfi_clustering(sdata_pp_bbknn, how="pbul", plot=False)
panfetal_utils.pfi_clustering(sdata_pp_bbknn, how="l", res=1.0, plot=False)
panfetal_utils.pfi_clustering(sdata_pp_bbknn, how="l", res=1.5, plot=False)

## Save output
adata_out_name = save_path + "{}.{}.batchCorrected.h5ad".format(suffix, s)
sdata_pp_bbknn.write_h5ad(adata_out_name)
