### READ OUTPUT FROM STARSOLO/CELLBENDER PIPELINE ###
import os, sys
import numpy as np
import scanpy as sc
import pandas as pd
import scrublet as scr
from multiprocessing import Pool
import tables
import numpy as np
import scipy.sparse as sp
from typing import Dict

## Reading cellbender outputs can have problem with 3GEX 
# Copy-pasting solution from https://github.com/broadinstitute/CellBender/issues/57#issuecomment-717370745
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.

    """

    try:
        import anndata
    except ImportError:
        raise ImportError('The anndata package must be installed to use the '
                          'function anndata_from_h5()')

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the count matrix.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)})
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata

## Function to make the output of anndata_from_h5 the same as sc.read_h5_10x
def uniform_output(adata):
    del adata.obs 
    del adata.uns
    del adata.obsm
    adata.var = adata.var[["id"]]
    adata.var.columns = ["gene_ids"]
    return adata

def read_cellbender_files(filename, raw_file_path, outdir,
               min_n_count = 2000, min_n_gene = 500, max_n_gene = 7000):
    ## Don't run if output exists already
    if "{s}_filtered.h5ad".format(s=filename) in os.listdir(outdir):
        print("Outfile exists already")
        return()
    else:
        ## Two-ways of loading cellbender outs
        path = '%s/%s/'%(raw_file_path,filename)
        try:
            adata = sc.read_10x_h5("{r}/{s}/{s}_filtered.h5".format(r=raw_file_path, s=filename))
        except:
            adata = anndata_from_h5("{r}/{s}/{s}_filtered.h5".format(r=raw_file_path, s=filename))
            adata = uniform_output(adata)

        adata.obs_names = [filename+"-"+x.strip("-1") for x in adata.obs_names]
        adata.var["GeneName"] = adata.var_names
        adata.var.columns = ["GeneID", "GeneName"]

        # caculate n_counts / n_genes per cell
        adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
        adata.obs['n_genes'] = np.sum(adata.X>0,axis=1)

        # filter cells
        print("Filtering cells...")
        clist = []
        clist.append(np.array(adata.obs['n_counts'] > min_n_count))
        clist.append(np.array(adata.obs['n_genes'] > min_n_gene))
        clist.append(np.array(adata.obs['n_genes'] < max_n_gene))

        c = np.column_stack(clist).all(axis=1)
        adata = adata[c].copy()

        adata = adata[:,np.argsort(adata.var.GeneID)]
        adata.obs['file'] = filename

        mito_genes = adata.var_names.str.startswith('MT-')
        adata.obs['mito'] = (np.sum(adata.X[:, mito_genes],axis=1).A1) / (np.sum(adata.X,axis=1).A1)

        print("Computing doublets...")
        scrub = scr.Scrublet(adata.X)
        if adata.shape[0] < 30:
            doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False, n_prin_comps=adata.shape[0] - 1)
        else:
            doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets
        sc.write('/%s/%s_filtered'%(outdir, filename),adata)
        
## dir of cellbender outputs
cellbender_dir = "/lustre/scratch117/cellgen/team205/sharedData/ly5/cellbender/"
## dir to save h5ad
outdir = '/nfs/team205/ed6/data/Fetal_immune/cellbender_raw/'
## Pick sample names from metadata table
metadata = pd.read_csv("/home/jovyan/mount/gdrive/Pan_fetal/annotations/manifest_clean_120121.csv", index_col=0)
sample_ls = metadata["Sample.lanes"][metadata["Sample.lanes"].isin(os.listdir(cellbender_dir))].tolist()

with Pool(10) as p:
    p.starmap(read_cellbender_files, [(s, cellbender_dir, outdir) for s in sample_ls])
