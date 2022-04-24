### Update scVI models ### 
## to make them compatible with latest scvi-tools version
import sys,os
import scvi
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy

import torch
device = torch.device("cuda")

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("m", help="model directory name")
parser.add_argument("spl", help="name of data split")
parser.add_argument("--model_dir", 
                    default='/home/jupyter/mount/gdrive/Pan_fetal/data4gpu_node/')
args = parser.parse_args()

m = args.m
spl = args.spl
model_dir = args.model_dir

adata = sc.read_h5ad(model_dir + f'PAN.A01.v01.entire_data_raw_count.20210429.{spl}.h5ad')
adata.var_names = adata.var['GeneID'].values.copy()

model_vars = pd.read_csv(model_dir + m + '/var_names.csv', header=None)[0]

adata_model = adata[:,model_vars].copy()

new_obs = pd.read_csv(model_dir + "PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.csv", index_col=0)
adata_model.obs = new_obs.loc[adata_model.obs_names]
adata_model.obs["bbk"] = adata_model.obs["method"] + adata_model.obs["donor"]

vae = scvi.model.SCVI.load(model_dir + m, adata=adata_model)

## Save updated model
vae.save(model_dir + m + "_scvi{v}".format(v=''.join(scvi.__version__.split("."))))