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
import splitting_rules

import argparse
parser = argparse.ArgumentParser()
# parser.add_argument("split_adata", help="AnnData object for PFI atlas subset")
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma...)")
parser.add_argument("--no_labels", dest="no_labels", action='store_true',
                    help="Avoids using labels for ridge regression")
parser.add_argument("--no_organ", dest="no_organ", action='store_true',
                    help="Avoids using organ covariate for ridge regression")
args = parser.parse_args()

## Load split data
save_path = '/nfs/team205/ed6/data/Fetal_immune/'
suffix = 'PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118'
s = args.split_name
no_labels = args.no_labels
no_organ = args.no_organ

adata_name = save_path + "{}.{}.h5ad".format(suffix, s)
sdata = sc.read_h5ad(adata_name)

# ## Add updated obs
new_obs = pd.read_csv("/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.full_obs.csv", index_col=0)
new_obs = new_obs.loc[sdata.obs_names]
new_obs["uniform_label_propagated"] = sdata.obs["uniform_label_propagated"]
sdata.obs = new_obs

## Preprocess (feature selection, scale)
sdata_pp = panfetal_utils.pfi_preprocess(sdata, how="", genesets2remove=[genes.cc_genes, genes.IG_genes, genes.TCR_genes])

## Ridge regression
start=time.time()
if no_organ:
    batch_keys=['method','donor']
else:
    batch_keys=['method','donor','organ']
if no_labels:
    conf_key = []
else:
    conf_key = ['uniform_label_propagated']
panfetal_utils.ridge_regression(sdata_pp, batch_key=batch_keys,confounder_key=conf_key)
v2_time = time.time() - start
print("Ridge regression runtime: ", str(v2_time))

sdata_pp.X = sdata_pp.layers['X_remain']
del sdata_pp.layers['X_remain']
del sdata_pp.layers['X_explained']

## BBKNN and clustering 
# Define "batch"
if no_organ:
    sdata_pp.obs['bbk'] = [x+y+z for x,y,z in zip(sdata_pp.obs['donor'],sdata_pp.obs['method'])]
else:
    sdata_pp.obs['bbk'] = [x+y+z for x,y,z in zip(sdata_pp.obs['donor'],sdata_pp.obs['method'],sdata_pp.obs['organ'])]

# Remove samples with less than 20 cells in split
sdata_pp_bbknn = sdata_pp[sdata_pp.obs['bbk'].isin(sdata_pp.obs["bbk"].value_counts().index[sdata_pp.obs["bbk"].value_counts() > 20])]
panfetal_utils.pfi_clustering(sdata_pp_bbknn, how="pbu", plot=False)
# panfetal_utils.pfi_clustering(sdata_pp_bbknn, how="l", res=1.0, plot=False)
panfetal_utils.pfi_clustering(sdata_pp_bbknn, how="l", res=1.5, plot=False)

## Save output
adata_out_name = save_path + "{}.{}.batchCorrected.h5ad".format(suffix, s)
if no_organ:
    adata_out_name = save_path + "{}.{}.batchCorrected.no_organ.h5ad".format(suffix, s)
sdata_pp_bbknn.write_h5ad(adata_out_name)

# ## Save output
# adata_out_name = save_path + "{}.{}.batchCorrected.h5ad".format(suffix, s)
# if no_organ:
#     adata_out_name = save_path + "{}.{}.batchCorrected.no_organ.h5ad".format(suffix, s)
# sdata_pp_bbknn.write_h5ad(adata_out_name)

## Assign cells to future split
if s in splitting_rules.splitting_scheme.keys():
    ## Assign updated labels to unlabelled cells based on KNN graph
    # panfetal_utils._propagate_labels(sdata_pp_bbknn, anno_col="uniform_label")
    panfetal_utils._propagate_labels(sdata_pp_bbknn, anno_col="uniform_label_lvl0")
    panfetal_utils._propagate_labels(sdata_pp_bbknn, anno_col="uniform_label_expanded_merged")

    # Add predicted lvl0 label based on most frequent cell type
    cl_counts = sdata_pp_bbknn.obs.reset_index()[['uniform_label_expanded_merged', 'leiden_150','index']] \
        .dropna() \
        .groupby(['leiden_150', 'uniform_label_expanded_merged']) \
        .count().fillna(0).reset_index() \
        .pivot(columns=['leiden_150'], index=['uniform_label_expanded_merged'])

    cl_frac = (cl_counts/cl_counts.sum(0))
    cl_frac.columns = ['index_' + str(x[1]) for x in cl_frac.columns]
    max_cl = cl_frac.max()
    top_3_ls = []
    for cl in cl_frac.columns:
        top_3 = cl_frac[cl_frac.index!="nan"].nlargest(1, cl)[cl].index[0]
        top_3_ls.append(top_3)

    pred_labels_df = pd.DataFrame(top_3_ls, columns=['leiden_150_pred_label'], index=[int(x.split("index_")[1]) for x in cl_frac.columns])

    pred_labels_df.index = pred_labels_df.index.astype('str')

    sdata_pp_bbknn.obs['leiden_150_pred_label'] = pred_labels_df.loc[sdata_pp_bbknn.obs['leiden_150']]['leiden_150_pred_label'].values

    ## Assign labels to new splits, according to name of starting split
    pfi_splitting_df = sdata_pp_bbknn.obs[['leiden_150_pred_label']].drop_duplicates()
    for sp in splitting_rules.splitting_scheme[s]:
        pfi_splitting_df["split_" + sp] = [sp if x in splitting_rules.splitting_labels[sp] else np.nan for x in pfi_splitting_df['leiden_150_pred_label']]


    ## Merge to anndata
    sdata_pp_bbknn.obs = sdata_pp_bbknn.obs[[x for x in sdata_pp_bbknn.obs.columns if "split" not in x]]
    obs_names = sdata_pp_bbknn.obs_names
    sdata_pp_bbknn.obs = sdata_pp_bbknn.obs.merge(pfi_splitting_df, on=['leiden_150_pred_label', 'leiden_150_pred_label'], how='left', indicator=False)
    sdata_pp_bbknn.obs_names = obs_names

    # Save
    split_df = sdata_pp_bbknn.obs[[x for x in sdata_pp_bbknn.obs.columns if "split_" in x]]
    split_df["uniform_label_propagated"] = sdata_pp_bbknn.obs["uniform_label_expanded_merged_propagated"]
    split_df.to_csv(save_path + "splits_info/{}.split_info.{}.csv".format(suffix,s))