### Build and store PAGA graphs and ForceAtlas2 visualizations 4 trajectory inference analysis ##
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.sparse
import anndata
from bbknn import bbknn

## import utils
cwd = '../../utils/'
sys.path.append(cwd)

import genes
import panfetal_utils

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma...)")
args = parser.parse_args()

## Directory to save figures
figdir = "/home/jovyan/mount/gdrive/Pan_fetal/Updates_and_presentations/figures/PAGA_graphs/"
if os.path.exists(figdir):
    sc.settings.figdir = figdir
else:
    os.mkdir(figdir)
    sc.settings.figdir = figdir
  
## Dir to save intermediate outputs
writedir = "/nfs/team205/ed6/data/Fetal_immune/PAGA_FA_outs/"
if not os.path.exists(writedir):
    os.mkdir(writedir)

def run_PAGA_FA(adata, key, split_name, anno_col='anno_lvl_2'):
    ## PAGA by annotated cell types
    sc.tl.paga(adata, groups = anno_col, neighbors_key = key)
    plt.rcParams["figure.figsize"] = [8,8]
    sc.pl.paga(adata, color='age', threshold=0.3, layout='fa', frameon=False, fontsize=7, node_size_scale=0.5,
               max_edge_width=0.6, colorbar=False, save="_" + key + "_anno.pdf")
    ## FA graph layout
    sc.tl.draw_graph(adata, init_pos='paga', maxiter=100)
    draw_graph_fa_file = "X_draw_graph_fa_paga.{s}.{k}.{a}.npy".format(s=split_name, k=key, a=anno_col)
    np.save(writedir + draw_graph_fa_file, adata.obsm['X_draw_graph_fa'])
    sc.pl.draw_graph(
        adata,
        color=[anno_col, 'age'],
        legend_loc='on data', legend_fontsize=10, save="PAGA_FA.{s}.{k}.{a}.pdf".format(s=split_name, k=key, a=anno_col))
    ## Get annotation centroid positions
    scatter_array = adata.obsm['X_draw_graph_fa']
    scatter_df = pd.DataFrame(scatter_array, columns=["x", "y"])
    scatter_df["color_source"] = adata.obs[anno_col].values
    all_pos = (
                scatter_df
                .groupby('color_source', observed=True)
                .median()
            )
    all_pos = all_pos.loc[adata.obs[anno_col].cat.categories]
    sc.pl.paga(adata, color='age', threshold=0.7, layout='fa', frameon=False, fontsize=7, node_size_scale=1,
               max_edge_width=0.9, colorbar=False, pos=all_pos.values,
              save="PAGA_FA.{s}.{k}.{a}.pdf".format(s=split_name, k=key, a=anno_col))

## Load anndata object + annotations
adata = panfetal_utils._load_split_and_annotation(args.split_name)
## Remove cells missing annotation
adata = adata[adata.obs["anno_lvl_2"] !='nan']

## Make trimmed BBKNN graph (k=50)
adata.uns["bbknn"] = adata.uns["neighbors"]
adata.obsp["bbknn_connectivities"] = adata.obsp["connectivities"]
adata.obsp["bbknn_distances"] = adata.obsp["distances"]

bbknn(adata, batch_key = "bbk", n_pcs=30, approx=True, trim=50)
adata.uns["bbknn50"] = adata.uns["neighbors"]
adata.obsp["bbknn50_connectivities"] = adata.obsp["connectivities"]
adata.obsp["bbknn50_distances"] = adata.obsp["distances"]

adata.write_h5ad(writedir + args.split_name + "_bbknn50.h5ad")

## Run PAGA and draw graph
run_PAGA_FA(adata, key="bbknn50", split_name=args.split_name)