### Merging reference and query and computing embedding with RAPIDS ###
## to run on GPU node
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

    ## Read reference
    ref_adata = sc.read_h5ad(ref_data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.{s}.h5ad'.format(t=timestamp, s=split))
    ref_embedding = np.load(ref_data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.{s}.scVI_out.npy'.format(t=timestamp, s=split))
    ref_adata.obsm["X_scvi"] = ref_embedding
    ref_adata.var_names = ref_adata.var.GeneID

    concat_adata = ref_adata.concatenate(query_adata_mapped, batch_key="dataset", batch_categories=["reference", "query"])
    return(concat_adata)

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

## Compute UMAP
print("Running KNN search...\n")
sc.pp.neighbors(merged_adata, n_neighbors=30, use_rep="X_scvi")
print("Running UMAP...\n")
sc.tl.umap(merged_adata, min_dist = 0.01, spread = 2)

## Save 
print("Saving file...\n")
merged_adata.write_h5ad(query_h5ad_file.split(".h5ad")[0] + ".withReference.h5ad")