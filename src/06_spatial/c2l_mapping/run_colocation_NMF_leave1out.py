### NMF factors cross-validation ###
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata

from cell2location.models.downstream import CoLocatedGroupsSklearnNMF

def _train_nmf(X_location_all):
    ## Train model
    nmf_model = CoLocatedGroupsSklearnNMF(n_fact=n_fact, X_data=X_location_all, n_iter=20000)
    nmf_model.fit()
    nmf_model.sample2df(node_name="nUMI_factors", ct_node_name="cell_type_factors")

    ## Get cell type proportions
    X_cts = nmf_model.cell_type_fractions.copy()
    X_cts.index = X_location_all.columns.str.strip("q95cell_abundance_w_sf_")
    X_cts.columns = [x[1] for x in X_cts.columns.str.split("mean_cell_type_factors")]

    ## Get location proportions
    X_locations = nmf_model.location_factors_df.copy()
    X_locations.index = X_location_all.index
    X_locations.columns = [x[1] for x in X_locations.columns.str.split("mean_nUMI_factors")]
    return(X_locations, X_cts)


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file", type=str)
args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad_file)

org = adata.obs['Organ'].unique()[0]
outdir = '/nfs/team205/ed6/data/Fetal_immune/spatial_analysis_outs/V2/'
n_fact = 10
sample_col = 'sample'

## Take out one sample
for s in adata.obs[sample_col].unique().tolist() + ['none']:
    X_location_all = adata[adata.obs[sample_col] != s].obsm['q95_cell_abundance_w_sf']
    X_locations, X_cts = _train_nmf(X_location_all)
    X_cts.to_csv(outdir + f"{org}_leaveout_{s}.celltype_fractions.csv")
    X_locations.to_csv(outdir + f"{org}_leaveout_{s}.locations.csv")