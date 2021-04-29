import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import anndata

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("scvi_embedding_file", help="path to .npy file of scVI embedding")
# parser.add_argument("h5ad_file", help="path to .h5ad file of anndata object")
args = parser.parse_args()

def embed_and_cluster_scvi(adata, emb_file):
    X_scVI_emb = np.load(emb_file)
    adata.obsm["X_scvi"] = X_scVI_emb
    print("Computing KNN graph...")
    sc.pp.neighbors(adata, use_rep = "X_scvi", n_neighbors = 30, key_added="scvi")
#     ## UMAP
#     print("Computing UMAP...")
#     sc.tl.umap(adata, min_dist = 0.01, spread = 2, neighbors_key="scvi")
#     np.save(emb_file.rstrip(".npy") + ".UMAP.npy", adata.obsm["X_umap"])
    ## Clustering 
    print("Clustering...")
    sc.tl.leiden(adata, resolution=1.5, key_added='leiden_150', n_iterations=5, neighbors_key="scvi")
    adata.obs[["leiden_150"]].to_csv(emb_file.rstrip(".npy") + ".clustering.csv")
    
# h5ad_file = args.h5ad_file
h5ad_file = '/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118.h5ad'

emb_file = args.scvi_embedding_file
# emb_file = "/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/PAN.A01.v01.entire_data_raw_count.wGut.scVI_out.5000HVGS.removeCC.keepTCRBCR.10ldims.npy"

## Read adata
print("Loading data...")
adata = sc.read_h5ad(h5ad_file)

## Run
embed_and_cluster_scvi(adata, emb_file)