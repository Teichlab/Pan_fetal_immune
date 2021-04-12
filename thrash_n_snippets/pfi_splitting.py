import os,sys
import numpy as np 
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata
import scipy
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix, find
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

from time import time
from datetime import datetime
import seaborn as sns

## import utils
cwd = '../utils/'
sys.path.append(cwd)

import genes
import panfetal_utils

import argparse
parser = argparse.ArgumentParser()
# parser.add_argument("split_adata", help="AnnData object for PFI atlas subset")
parser.add_argument("--no_labels", dest="no_labels", action='store_true',
                    help="Avoids using labels for ridge regression")
parser.add_argument("--no_organ", dest="no_organ", action='store_true',
                    help="Avoids using organ covariate for ridge regression")

args = parser.parse_args()

## Set directory to save figures
figdir = "/home/jovyan/mount/gdrive/Pan_fetal/Updates_and_presentations/figures/clustering_and_splitting/"
if os.path.exists(figdir):
    sc.settings.figdir = "/home/jovyan/mount/gdrive/Pan_fetal/Updates_and_presentations/figures/clustering_and_splitting/"
else:
    os.mkdir(figdir)
    sc.settings.figdir = "/home/jovyan/mount/gdrive/Pan_fetal/Updates_and_presentations/figures/clustering_and_splitting/"
    

## Read data
adata = sc.read_h5ad("/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118.h5ad")

## Add clustering results
cl_obs = pd.read_csv("/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118.clustering.obs.csv", index_col=0)
cl_obs.index = adata.obs_names
adata.obs = pd.concat([adata.obs, cl_obs.loc[adata.obs_names][["leiden_100", "leiden_150"]]],axis=1)

## Some fixes to labels
## wrong labelling of Enterocytes and a typo
obs_df = adata.obs
obs_df["uniform_label_expanded_merged"] = obs_df["uniform_label_expanded_merged"].astype("str")
obs_df.loc[adata.obs_names[np.where(adata.obs.uniform_label=="ENTEROCYTE")], "uniform_label_expanded_merged"] = "ENTEROCYTE"
obs_df.loc[adata.obs_names[np.where(adata.obs.uniform_label=="LMPP")],"uniform_label_expanded_merged"] = "LYMPHOID PROGENITOR"
adata.obs = obs_df

## Propagate labels to cells with missing labels (by KNN mean)
panfetal_utils._propagate_labels(adata, anno_col="uniform_label")
panfetal_utils._propagate_labels(adata, anno_col="uniform_label_lvl0")
panfetal_utils._propagate_labels(adata, anno_col="uniform_label_expanded_merged")

sc.pl.umap(adata, color="organ", save='_organ.png')
sc.pl.umap(adata, color="leiden_150_pred_label_expanded", legend_loc="on data", save='_leiden_label.png')
sc.pl.umap(adata, color=["uniform_label_propagated"],palette=sc.plotting.palettes.default_102, save='_propagated_labels.png')

## Add putative labels based on the most abundant cells in clusters
adata.obs['leiden_150'] = [str(x) for x in adata.obs['leiden_150']]

## Add predicted lvl0 label based on most frequent cell type
cl_counts = adata.obs.reset_index()[['uniform_label_lvl0_propagated', 'leiden_150','index']] \
    .dropna() \
    .groupby(['leiden_150', 'uniform_label_lvl0_propagated']) \
    .count().fillna(0).reset_index() \
    .pivot(columns=['leiden_150'], index=['uniform_label_lvl0_propagated'])

cl_frac = (cl_counts/cl_counts.sum(0))
cl_frac.columns = ['index_' + str(x[1]) for x in cl_frac.columns]
max_cl = cl_frac.max()
top_3_ls = []
for cl in cl_frac.columns:
    top_3 = cl_frac[cl_frac.index!="nan"].nlargest(1, cl)[cl].index[0]
    top_3_ls.append(top_3)

pred_labels_df = pd.DataFrame(top_3_ls, columns=['leiden_150_pred_label'], index=[int(x.split("index_")[1]) for x in cl_frac.columns])

pred_labels_df.index = pred_labels_df.index.astype('str')

adata.obs['leiden_150_pred_label'] = pred_labels_df.loc[adata.obs['leiden_150']]['leiden_150_pred_label'].values

## Add predicted uniform label based on most frequent cell type
cl_counts = adata.obs.reset_index()[["uniform_label_expanded_merged_propagated", 'leiden_150','index']] \
    .dropna() \
    .groupby(['leiden_150', "uniform_label_expanded_merged_propagated"]) \
    .count().fillna(0).reset_index() \
    .pivot(columns=['leiden_150'], index=["uniform_label_expanded_merged_propagated"])

cl_frac = (cl_counts/cl_counts.sum(0))
cl_frac.columns = ['index_' + str(x[1]) for x in cl_frac.columns]
max_cl = cl_frac.max()
top_3_ls = []
for cl in cl_frac.columns:
    top_3 = cl_frac[cl_frac.index!="nan"].nlargest(1, cl)[cl].index[0]
    top_3_ls.append(top_3)

pred_labels_df = pd.DataFrame(top_3_ls, columns=['leiden_150_pred_label'], index=[int(x.split("index_")[1]) for x in cl_frac.columns])

pred_labels_df.index = pred_labels_df.index.astype('str')
# pred_labels_df.loc[adata.obs['leiden_100']]['leiden_100_pred_label']

adata.obs['leiden_150_pred_label_expanded'] = pred_labels_df.loc[adata.obs['leiden_150']]['leiden_150_pred_label'].values