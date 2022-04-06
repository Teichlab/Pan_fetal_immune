import sys,os
import scvi
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy

cwd = '.'
sys.path.append(cwd)
import map_query_utils

import torch
device = torch.device("cuda")

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file", help="path to query anndata file")
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma...)")
parser.add_argument("batch_name", 
                    help='''
                    how do you want to label the new data you are adding? Used only if query_adata.obs['bbk'] is None
                    ''')
parser.add_argument("--indir", 
                    default="/home/jupyter/mount/gdrive/Pan_fetal/data4gpu_node/",
                    help="folder containing scvi outputs")
parser.add_argument("--model_dir", 
                    default=None,
                    help="folder containing scvi model (Default: None, uses scvi_{split}_model/)")
parser.add_argument("--query_gene_id_col", 
                    default="gene_ids",
                    help="folder containing scvi outputs")
parser.add_argument("--w_decay", 
                    default=0.0, type=float,
                    help="Weight decay parameter (default: 0.0, no parameter updating)")
args = parser.parse_args()

h5ad_file = args.h5ad_file
split = args.split_name
batch_name = args.batch_name
scvi_outs_dir = args.indir
model_dir = args.model_dir
w_decay = args.w_decay

def _verify_counts(adata):
    return(all([not (i%1) for i in adata.X[0,:].toarray()[0]]))

def _map_query_to_panfetal_v2(
    query_adata,
    split,
    scvi_outs_dir = "/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/",
    model_dir = None,
    timestamp = "20210429",
    w_decay = 0.0
    ):
    '''
    Use scArches method to do online update of scVI model with query dataset
    
    Params:
    -------
    - query_adata: anndata object of query dataset to map
    - split: name of split to use (to chose the model)
    - scvi_outs_dir: directory storing the scVI training outputs
    - model_dir: which model directory in `scvi_outs` to use? (Default: None, uses `scvi_{split}_model/`)
    
    Outputs:
    - anndata with new embedding in obsm ("X_scVI_project")
    - vae: model
    '''
    if model_dir is None:
        model_dir='scvi_' + split + '_model/'
    ## Check that the feature used for reference training match the var_names in query data
    var_names_model = torch.load(scvi_outs_dir + model_dir + 'model.pt')['var_names']

    vars_overlap = any(pd.Series(var_names_model).isin(query_adata.var_names))

    ## Extract gene ids if not
    if not vars_overlap:
        raise ValueError("Vars don't overlap -- convert to geneIDs")

    # ## Add batch column
    # if "bbk" not in query_adata.obs.columns:
    #     print("considering all query cells from the same technical batch")
    #     query_adata.obs["bbk"] = batch_name

    ## Zero-fill missing genes in query
    is_in_query_var = pd.Series(var_names_model).isin(query_adata.var_names)
    n_genes = len(var_names_model[~is_in_query_var])
    if n_genes > 0:
        empty_X = np.zeros(shape=[ query_adata.n_obs, n_genes])
        empty_query_adata = anndata.AnnData(X=empty_X, obs=query_adata.obs)
        empty_query_adata.var_names = var_names_model[~is_in_query_var]
        empty_query_adata.var_names.names = ["index"]
        query_adata_filled = anndata.concat([query_adata, empty_query_adata], axis=1)
        query_adata_filled = query_adata_filled[:,var_names_model].copy()
        query_adata_filled.obs = query_adata.obs.copy()
    else:
        query_adata_filled = query_adata.copy()

    ## Load new model with the query data
    vae_q = scvi.model.SCVI.load_query_data(
        query_adata_filled,
        scvi_outs_dir + model_dir,
        inplace_subset_query_vars=True
    )

    ## Train
    vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=w_decay))
    query_adata_filled.obsm["X_scvi"] = vae_q.get_latent_representation()
    return(query_adata_filled, vae_q)


query_adata = sc.read_h5ad(h5ad_file)

## Check that varnames are gene ids
if not query_adata.var_names.str.startswith("ENS").all():
    raise ValueError('`query_adata.var_names` are not Ensembl geneIDs. Please convert')

## Check that adata.X contains counts
if not _verify_counts(query_adata):
    raise ValueError('`query_adata.X` does not contain raw counts')

query_adata_mapped, _ = _map_query_to_panfetal_v2(query_adata, split=split,  
                                            scvi_outs_dir = scvi_outs_dir,
                                            model_dir=model_dir,
                                                  w_decay=w_decay
                                                 )

if w_decay==0.0:
    out_file = '{p}.mapped2{s}.h5ad'.format(p=h5ad_file.strip(".h5ad"), s=split)
else:
    out_file = '{p}.mapped2{s}.nonzero_w_decay.h5ad'.format(p=h5ad_file.strip(".h5ad"), s=split)
query_adata_mapped.write_h5ad(out_file)
