### NMF factors cross-validation ###
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata

from cell2location.models.downstream import CoLocatedGroupsSklearnNMF

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("organ", type=str)
parser.add_argument("n_fact", type=int) 
parser.add_argument("--filter_cts",action='store_true',
                    help="Should only common cell types measured early and late be kept")
parser.add_argument("--outdir",
                    default='/nfs/team205/ed6/data/Fetal_immune/spatial_analysis_outs/V2/',
                    help="path to working directory")
args = parser.parse_args()

org = args.organ
outdir = args.outdir
n_fact = args.n_fact
filter_cts = args.filter_cts

adata_early = sc.read_h5ad(outdir + f"{org}_early.h5ad")
adata_late = sc.read_h5ad(outdir + f"{org}_late.h5ad")

## Merge in one obj for NMF analysis
X_location_early = adata_early.obsm['q95_cell_abundance_w_sf'].copy()
X_location_late = adata_late.obsm['q95_cell_abundance_w_sf'].copy()
X_location_all = pd.concat([X_location_early, X_location_late])
if filter_cts:
    ## Remove stage specific cell types
    X_location_all = pd.concat([X_location_early, X_location_late])
    X_location_all = X_location_all.loc[:,~X_location_all.isna().any()] 
else:
    X_location_all = X_location_all.fillna(0) ## Set missing cts to 0

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

adata_early.obsm['nmf'] = X_locations.loc[adata_early.obs_names]
adata_late.obsm['nmf'] = X_locations.loc[adata_late.obs_names]

## Save

if filter_cts:
    adata_early.write_h5ad(outdir + f"{org}_early.filter_cts.h5ad")
    adata_late.write_h5ad(outdir + f"{org}_late.filter_cts.h5ad")
    X_cts.to_csv(outdir + f"{org}_celltype_fractions.filter_cts.csv")
else:
    adata_early.write_h5ad(outdir + f"{org}_early.h5ad")
    adata_late.write_h5ad(outdir + f"{org}_late.h5ad")
    X_cts.to_csv(outdir + f"{org}_celltype_fractions.csv")

# ## Run on farm
# conda activate scvi-env
# echo "python run_colocation_NMF.py TH ${n}" | bsub -G teichlab -o logfile-%J.txt -M20000 -R "select[mem>20000] rusage[mem=20000]"
