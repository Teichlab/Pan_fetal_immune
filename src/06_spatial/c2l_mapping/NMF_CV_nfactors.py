### NMF factors cross-validation ###
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata
import pickle

from cell2location.models.downstream import CoLocatedGroupsSklearnNMF
from sklearn.metrics import mean_squared_error

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("organ", type=str)
parser.add_argument("fold_id", type=str) 
parser.add_argument("n_fact", type=int) 
parser.add_argument("--outdir",
                    default='/nfs/team205/ed6/data/Fetal_immune/spatial_analysis_outs/V2/',
                    help="path to working directory")
args = parser.parse_args()

def _split_for_cv(X_location_all, frac_cv = 0.05, seed=42):
    np.random.seed(seed)
    size_cv = np.round(X_location_all.shape[0]*frac_cv)
    cv_ixs = np.random.choice(np.arange(X_location_all.shape[0]), size=int(size_cv))
    cv_ixs.sort()

    cv_spots = X_location_all.index[cv_ixs]
    train_spots = X_location_all.index[~X_location_all.index.isin(cv_spots)]
    if not len(np.intersect1d(train_spots, cv_spots)) == 0:
        print('Something going wrong')

    X_train = X_location_all.loc[train_spots]
    X_cv = X_location_all.loc[cv_spots]
    return(X_train, X_cv)

def _cv_model(X_location_all, n_fact, frac_cv = 0.05, n_iter=10, seed=42):
    ## Split data
    X_train, X_cv = _split_for_cv(X_location_all, frac_cv=frac_cv, seed=seed)
    ## Fit model to train
    nmf_model = CoLocatedGroupsSklearnNMF(n_fact=n_fact, X_data=X_train, n_iter=n_iter)
    nmf_model.fit(n=1)

    ## Predict on CV dataset
    pred_factors = nmf_model.models['init_1'].transform(X_cv)

    ## Compute reconstruction error
    W = nmf_model.models['init_1'].components_.copy()
    recon_cv_X = pred_factors.dot(W)
    mse = mean_squared_error(recon_cv_X.flatten(), X_cv.values.flatten())
    return(mse)

org = args.organ
fold_id = args.fold_id
outdir = args.outdir
n_fact = args.n_fact

adata_early = sc.read_h5ad(outdir + f"{org}_early.h5ad")
adata_late = sc.read_h5ad(outdir + f"{org}_late.h5ad")

## Merge in one object for NMF analysis
X_location_early = adata_early.obsm['q95_cell_abundance_w_sf'].copy()
X_location_late = adata_late.obsm['q95_cell_abundance_w_sf'].copy()
X_location_all = pd.concat([X_location_early, X_location_late]).fillna(0) ## Set missing cts to 0

mse = _cv_model(X_location_all, n_fact, n_iter=15000, seed=int(fold_id))
with open(outdir + f"MSE_{org}_fold{fold_id}_n_fact{n_fact}.pkl", 'wb') as f:
    pickle.dump((org, fold_id, n_fact, mse), f)
    
# ## Run on farm
# conda activate scvi-env
# for f in $(seq 0 5); do
#     for n in $(seq 5 20); do
#         echo "python NMF_CV_nfactors.py TH ${f} ${n}" | bsub -G teichlab -o logfile-%J.txt -M20000 -R "select[mem>20000] rusage[mem=20000]"
#         done
#     done