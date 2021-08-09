### Merging reference and query and computing embedding ###

import os,sys
import numpy as np
import scanpy as sc
import anndata

def _merge_query_and_reference(
    query_h5ad_file,
    split,
    ref_data_dir = '/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/',
    timestamp = '20210429'):
    '''
    Read output from query mapping and reference data and merges them in one anndata
    
    Params:
    ------
    - query_h5ad_file: h5ad file storing mapped query anndata
    - split: data split ID
    - ref_data_dir: folder containing reference data
    - timestamp: dataset timestamp
    '''
    query_adata_mapped = sc.read_h5ad(query_h5ad_file)
    query_adata_full = sc.read_h5ad(query_h5ad_file.split(".mapped2")[0] + ".h5ad") ## To add genes that are not used in scVI
    query_adata_full.obsm["X_scvi"] = query_adata_mapped.obsm["X_scvi"].copy()
    query_adata_full.uns["_scvi"] = query_adata_mapped.uns["_scvi"].copy()
    
    ## Read reference
    ref_adata = sc.read_h5ad(ref_data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.{s}.h5ad'.format(t=timestamp, s=split))
    ref_embedding = np.load(ref_data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.{s}.scVI_out.npy'.format(t=timestamp, s=split))
    ref_adata.obsm["X_scvi"] = ref_embedding
    ref_adata.var_names = ref_adata.var.GeneID

    concat_adata = anndata.concat([ref_adata, query_adata_full], axis=0,
                                  label="dataset", keys=["reference", "query"],
                                  join="outer", merge="unique", uns_merge="unique")
    concat_adata.obs_names = concat_adata.obs_names + "-" + concat_adata.obs["dataset"].astype("str")
    return(concat_adata)

def _add_all_query_genes(merged_adata, query_adata_full):
    if not any(merged_adata.var_names.isin(query_adata_full.var_names)):
        raise ValueError("var_names don't match between query and merged AnnData")

    if not any(query_adata_full.obs_names.str.endswith("-query")):
        query_adata_full.obs_names = query_adata_full.obs_names + "-query"

    ## Do the merge
    full_merged_adata = anndata.concat([merged_adata, query_adata_full], axis=1, join="outer", merge="unique", uns_merge="unique")

    ## Check that the number of obs is right
    if not full_merged_adata.n_obs == merged_adata.n_obs:
        raise AssertionError("The number of obs doesn't match, something is wrong in your join")

    ## Check that you have more expression vals than before
    one_query_cell = query_adata_full.obs_names[10]
    if not len(full_merged_adata[one_query_cell].X.nonzero()[0]) > len(merged_adata[one_query_cell].X.nonzero()[0]):
        raise AssertionError("You have less or the same expression values for query cells than before. Are you sure that query_adata_full is the dataset before feature selection?")

    return(full_merged_adata)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("query_h5ad_file", help="path to query anndata file")
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma...)")
parser.add_argument("--ref_data_dir", 
                    default="/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/",
                    help="folder containing reference data")
parser.add_argument("--timestamp", 
                    default="20210429",
                    help="data time stamp")
args = parser.parse_args()

query_h5ad_file = args.query_h5ad_file
split = args.split_name
ref_data_dir = args.ref_data_dir
timestamp = args.timestamp

## Merge datasets
print("Merging reference and query...\n")
merged_adata = _merge_query_and_reference(query_h5ad_file, split, ref_data_dir=ref_data_dir)
query_adata_full = sc.read_h5ad(query_h5ad_file.split(".mapped2")[0] + ".h5ad") ## To add genes that are not used in scVI
# merged_adata = _add_all_query_genes(merged_adata, query_adata_full)

## Compute UMAP
print("Running KNN search...\n")
sc.pp.neighbors(merged_adata, n_neighbors=30, use_rep="X_scvi")
print("Running UMAP...\n")
sc.tl.umap(merged_adata, min_dist = 0.01, spread = 2)

## Save 
print("Saving file...\n")
merged_adata.write_h5ad(query_h5ad_file.split(".h5ad")[0] + ".withReference.h5ad")