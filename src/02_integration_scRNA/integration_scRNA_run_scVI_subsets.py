#!/usr/bin/env python

## Run scVI on dataset splits

import sys,os
import scvi
import anndata
import matplotlib
import seaborn as sns
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import scanpy as sc
import numpy.random as random

import torch
device = torch.device("cuda")

cwd = '../utils/'
sys.path.append(cwd)
import genes
import panfetal_utils

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("split_name", 
                    default="",
                    help="ID for data split (e.g. NKT, Progenitors, Stroma...) (default: no split, full atlas)")
args = parser.parse_args()
spl = args.split_name

def load_data_split(data_dir, timestamp, split):
    # Load estimated gene dispersions for HVG selection
    # Generated running `Pan_fetal_immune/utils/PFI_pp_4_HVG_stats.py`
    adata_lognorm_var = pd.read_csv(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.var.csv'.format(t=timestamp, s=split))

    ### Load count data
    adata_raw = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.{s}.h5ad'.format(t=timestamp, s=split))
    adata_raw.var_names_make_unique()

    # Load obs
    new_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.csv".format(t=timestamp), index_col=0)
    adata_raw.obs = new_obs.loc[adata_raw.obs_names]

    ## Load var
    adata_raw.var = adata_lognorm_var.copy()

    ## Add batch key
    adata_raw.obs["bbk"] = adata_raw.obs["method"] + adata_raw.obs["donor"]
    return(adata_raw)

def subset_top_hvgs(adata_lognorm, n_top_genes):
    dispersion_norm = adata_lognorm.var['dispersions_norm'].values.astype('float32')

    dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
    dispersion_norm[
                ::-1
            ].sort()  # interestingly, np.argpartition is slightly slower

    disp_cut_off = dispersion_norm[n_top_genes - 1]
    gene_subset = adata_lognorm.var['dispersions_norm'].values >= disp_cut_off
    return(adata_lognorm[:,gene_subset])

def prep_scVI(adata, 
              n_hvgs = 5000,
              remove_cc_genes = True,
              remove_tcr_bcr_genes = False
             ):
    ## Remove cell cycle genes
    if remove_cc_genes:
        adata = panfetal_utils.remove_geneset(adata,genes.cc_genes)

    ## Remove TCR/BCR genes
    if remove_tcr_bcr_genes:
        adata = panfetal_utils.remove_geneset(adata, genes.IG_genes)
        adata = panfetal_utils.remove_geneset(adata, genes.TCR_genes)
        
    ## HVG selection
    adata = subset_top_hvgs(adata, n_top_genes=n_hvgs)
    return(adata)

def train_scVI(adata, n_dims=20):
    adata = scvi.data.setup_anndata(adata, batch_key = "bbk", copy=True)
    arches_params = dict(
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
        )
    vae = scvi.model.SCVI(adata, n_latent=n_dims, **arches_params)
    vae.train(early_stopping=True,
        train_size=0.9,
        early_stopping_patience=45,
        max_epochs=400, 
        batch_size=1024, 
        limit_train_batches=20
       )
    return(vae)

def scvi_split(s, 
               data_dir = "/home/jupyter/mount/gdrive/Pan_fetal/data4gpu_node/",
               figdir = '/home/jupyter/mount/gdrive/Pan_fetal/Updates_and_presentations/',
               timestamp = "20210429"
              ):
    adata_raw = load_data_split(data_dir, timestamp, s)
    adata_raw = prep_scVI(adata_raw, n_hvgs=7500, remove_cc_genes=True, remove_tcr_bcr_genes=True)
    vae = train_scVI(adata_raw, n_dims=20)
    adata_raw.obsm["X_scVI"] = vae.get_latent_representation()
    ## Save embedding
    outname = "PAN.A01.v01.entire_data_raw_count.{t}.{s}.scVI_out.npy".format(t=timestamp, s=s)
    np.save(data_dir + outname, adata_raw.obsm["X_scVI"])
    ## Plot convergence
    sns.set_context("talk")
    plt.plot(vae.history["elbo_train"], label="Training");
    plt.plot(vae.history["elbo_validation"], label="Validation");
    plt.legend();
    plt.xlabel("epoch");
    plt.ylabel("ELBO");
    plt.savefig("{f}/scvi_training_elbo_{s}.pdf".format(s=s, f=figdir), bbox_inches="tight")
    # save the reference model
    model_dir = 'scvi_{s}_model/'.format(s=s)
    if not os.path.exists(data_dir + model_dir):
        os.mkdir(data_dir + model_dir)
    vae.save(data_dir + model_dir, overwrite=True)

scvi_split(spl)
