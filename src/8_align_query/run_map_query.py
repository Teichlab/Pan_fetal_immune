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
args = parser.parse_args()

h5ad_file = args.h5ad_file
split = args.split_name
batch_name = args.batch_name
scvi_outs_dir = args.indir
model_dir = args.model_dir

def _verify_counts(adata):
    return(all([not (i%1) for i in adata.X[0,:].toarray()[0]]))

query_adata = sc.read_h5ad(h5ad_file)

## Check that varnames are gene ids
if not query_adata.var_names.str.startswith("ENS").all():
    raise ValueError('`query_adata.var_names` are not Ensembl geneIDs. Please convert')
    
## Check that adata.X contains counts
if not _verify_counts(query_adata):
    raise ValueError('`query_adata.X` does not contain raw counts')
    
query_adata_mapped, _ = map_query_utils._map_query_to_panfetal(query_adata, split=split,  
                                            scvi_outs_dir = scvi_outs_dir,
                                            batch_name=batch_name, model_dir=model_dir)

out_file = '{p}.mapped2{s}.h5ad'.format(p=h5ad_file.strip(".h5ad"), s=split)
query_adata_mapped.write_h5ad(out_file)