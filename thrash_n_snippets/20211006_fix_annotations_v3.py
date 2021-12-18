#!/usr/bin/env python
import os,sys
import numpy as np 
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata
import scipy
import matplotlib.pyplot as plt
import re
import scvelo as scv
import seaborn as sns

data_dir = '/nfs/team205/ed6/data/Fetal_immune/'
timestamp = '20210429'
marker_genes_dir = "/home/jovyan/mount/gdrive/Pan_fetal/significant_genes/marker genes/"

figdir = "/home/jovyan/mount/gdrive/Pan_fetal/Updates_and_presentations/figures/annotationsV2/"
if os.path.exists(figdir):
    sc.settings.figdir = figdir
    scv.settings.figdir = figdir
else:
    os.mkdir(figdir)
    sc.settings.figdir = figdir
    scv.settings.figdir = figdir

def _load_split(split, min_n_cells = 50):
    spl_adata = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.embedding.h5ad'.format(t=timestamp, s=split))
    spl_adata = spl_adata[~spl_adata.obs_names.isin(mat_barcodes)]
    
    ## Save annotations for split
    spl_adata.obs["anno_lvl_2_final_clean"] = anno_obs.loc[spl_adata.obs_names]["anno_lvl_2_final_clean"].astype("str")

    ## Label as OTHER cell labels that appear in less than 50 cells
    small_labels = spl_adata.obs["anno_lvl_2_final_clean"].value_counts()[spl_adata.obs["anno_lvl_2_final_clean"].value_counts() < min_n_cells].index
    spl_adata.obs.loc[spl_adata.obs["anno_lvl_2_final_clean"].isin(small_labels), "anno_lvl_2_final_clean"] = "OTHER"
    spl_adata.obs.loc[spl_adata.obs["anno_lvl_2_final_clean"]=="LOW_Q_INCONSISTENT", "anno_lvl_2_final_clean"] = "OTHER"
    return(spl_adata)

def _plot_marker_dotplot(spl_adata, mark_split, 
                         labels_order=None, 
                         gene_labels_order=None, 
                         unique_genes=False,
                         **kwargs):
    markers_df = all_markers_df[all_markers_df["split"].isin(mark_split)]
    markers_df = markers_df[markers_df.gene.isin(spl_adata.var_names)]
    # markers_df = markers_df.sort_values('anno_lvl_2')
    if gene_labels_order is None:
        gene_labels_order = markers_df.anno_lvl_2.unique().tolist()

    markers_df= markers_df[markers_df.anno_lvl_2.isin(gene_labels_order)]
    markers_dict = {g: d['gene'].values.tolist() for g, d in markers_df.groupby('anno_lvl_2')}
    markers_dict = {k: markers_dict[k] for k in gene_labels_order}
    
    if labels_order is None:
        labels_order = gene_labels_order
    
    keep_labels = [x for x in labels_order if x in spl_adata.obs["anno_lvl_2_final_clean"].values]
    
    if unique_genes:
        markers_ls_2 = sum(markers_dict.values(), [])
        markers_ls_ixs = np.unique(markers_ls_2, return_index=True)[1]
        markers_dict = [markers_ls_2[index] for index in sorted(markers_ls_ixs)]
    
    sc.pl.dotplot(spl_adata[spl_adata.obs["anno_lvl_2_final_clean"].isin(keep_labels)], 
                      markers_dict, 
                      groupby="anno_lvl_2_final_clean", 
                      categories_order=keep_labels, **kwargs)
        
        
def _plot_split_embedding(spl_adata, split, **kwargs):
    sns.set_context("talk")
    scv.pl.umap(spl_adata, color="anno_lvl_2_final_clean", legend_loc="on data", 
                size=2, legend_fontoutline=3,
                title=split, **kwargs)


def _save_plot_markers_embedding(spl_adata, mark_split, ncols=4, gene_labels_order = None):
    markers_df = all_markers_df[all_markers_df["split"].isin(mark_split)]
    markers_df = markers_df[markers_df.gene.isin(spl_adata.var_names)]
    # markers_df = markers_df.sort_values('anno_lvl_2')
    if gene_labels_order is None:
        gene_labels_order = markers_df.anno_lvl_2.unique().tolist()
    emb_figdir = "embedding_markers_{s}/".format(s=''.join(mark_split))
    if not os.path.exists(figdir + emb_figdir):
        os.mkdir(figdir + emb_figdir)
    scv.settings.figdir = figdir + emb_figdir

    for label in markers_df.anno_lvl_2.unique():
        markers_label_df = markers_df[markers_df.anno_lvl_2==label]
        sns.set_context("talk")
        marker_genes = markers_label_df.gene.tolist()
        scv.pl.umap(spl_adata, color=marker_genes, size=5, fontsize=32, legend_fontsize=32, ncols=ncols,
                    show=False,
                   save="embedding_markers_anno_lvl_2_{s}.{a}.png".format(s="".join(mark_split), a="_".join(label.split("/"))))
    scv.settings.figdir = figdir

# get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# # In[9]:


# get_ipython().run_cell_magic('R', '', 'library(tidyverse)\nlibrary(ggplot2)')


# ## Build annotation file
# Load annotation file curated by Chenqu/Isaac/Laura
anno_obs_new = pd.read_csv('../../metadata/manual_annotation/IG_anno_lvl_2_final_clean_051121.csv', index_col=0)

# Load old annotation file
anno_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.clean.csv".format(t=timestamp), index_col=0)
# anno_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.clean.V2.csv".format(t=timestamp), index_col=0)

# ## Save 4 legacy
# anno_obs.to_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.clean.V2.csv".format(t=timestamp))


# In[119]:


anno_obs.anno_lvl_2_final_clean = anno_obs_new.anno_lvl_2_final_clean.copy()


# In[120]:


anno_obs = anno_obs[['n_counts', 'n_genes', 'file', 'mito', 'doublet_scores',
       'predicted_doublets', 'name', 'uniform_label',
       'uniform_label_expanded_merged', 'uniform_label_lvl0', 'organ',
       'Sample.lanes', 'Sort_id', 'age', 'method', 'donor', 'sex',
       'Processing_method', 'AnnatomicalPart', 'Sample', 'anno_lvl_2_final_clean']]


# In[14]:


all_labels = anno_obs.anno_lvl_2_final_clean.unique()

## clean
clean = {
    ## Myeloid labels
    'LC':'LANGERHANS_CELLS',
    "MONOCYTE_III_CCR2":'MONOCYTE_II_CCR2',
    "MONOCYTE_II_IL1B":"MONOCYTE_III_IL1B",
    'MACROPHAGE_I':"MACROPHAGE_LYVE1_HIGH",
    'MACROPHAGE_II':"MACROPHAGE_IRON_RECYCLING",
    'MACROPHAGE_III':"MACROPHAGE_PROLIFERATING",
    "MACROPHAGE_IV":"MACROPHAGE_MHCII_HIGH",
    'MACROPHAGE_V':"MACROPHAGE_KUPFFER_LIKE",
    'MACROPHAGE_VIII_TREM2':'MACROPHAGE_TREM2',
    'EOSINOPHILBASOPHIL':'EOSINOPHIL_BASOPHIL',
    'MACROPHAGE_VII_ERY':'MACROPHAGE_ERY',
    'MACROPHAGE_VI_PERI':'MACROPHAGE_PERI',
    ## Lymphoid labels
    'TH17': 'TYPE_3_INNATE_T',
    'NK_T': 'TYPE_1_INNATE_T',
    'PRE_PRO_B_CELL':'PRE_PRO_B',
    'DN(P)_T_CELL':'DN(P)_T',
    'DN(early)_T_CELL':'DN(early)_T',
    'ELP':"PRE_PRO_B"
    }

clean_anno = [clean[x] if x in clean.keys() else x for x in anno_obs.anno_lvl_2_final_clean]
clean_anno = ["_".join(x.split("/")) if "/" in x else x for x in clean_anno]
clean_anno = [re.split(" CELLS?$", x)[0] for x in clean_anno]
clean_anno = ["_".join(x.split(" ")) for x in clean_anno]
clean_anno = pd.Series(clean_anno)

anno_obs.anno_lvl_2_final_clean = clean_anno.values.copy()


# #### Update groupings

# In[15]:


## Read annotation groupings
import json
with open('../../metadata/anno_groups.json', 'r') as json_file:
    anno_groups_dict = json.load(json_file)

anno_groups_dict_rev = {x:g for g,a in anno_groups_dict.items() for x in a}


# In[16]:


all_labels = anno_obs.anno_lvl_2_final_clean.unique()

missing_group = [x for x in all_labels if x not in anno_groups_dict_rev.keys()]
for k in missing_group:
    print(k)
    if "FIBROBLAST" in k or "ENDOTHEL" in k or "HEPATOCYTE" in k or "VSMC" in k or "EPITHELIUM" in k or "ENTEROENDOCRINE" in k:
        anno_groups_dict_rev[k] = "STROMA"
    if k=='LMPP_MLP':
        anno_groups_dict_rev[k] = "PROGENITORS"
    if k.endswith("_T"):
        anno_groups_dict_rev[k] = "NK/T CELLS"
    else:
        anno_groups_dict_rev[k] = "MYELOID"

## Remove old labels        
anno_groups_dict_rev = {k:v for k,v in anno_groups_dict_rev.items() if k in all_labels}

ks = list(set([x for x in anno_groups_dict_rev.values()]))
anno_groups_dict = {}
for k1 in ks:
    anno_groups_dict[k1] = [k for k,v in anno_groups_dict_rev.items() if anno_groups_dict_rev[k]==k1]


# In[17]:


with open('../../metadata/anno_groups.json', 'w') as outfile:
    json.dump({k:list(v) for k,v in anno_groups_dict.items()}, outfile)


# ## Compare published VS new annotations

# In[58]:


reannotated_df = anno_obs[~anno_obs['uniform_label'].isna()]
annotations_conf_mat = sc.metrics.confusion_matrix("uniform_label", "anno_lvl_2_final_clean", reannotated_df, normalize=False)
anno_tot=sc.metrics.confusion_matrix("uniform_label", "anno_lvl_2_final_clean", anno_obs[~anno_obs['uniform_label'].isna()], normalize=False).sum(0)


# In[46]:


# ## Group labels for plotting
# progenitors_labels = [x for x in reannotated_df.anno_lvl_0_final.unique() if x not in anno_vocabulary["anno_lvl_0"].unique() and x != "nan" and "T CELL" not in x]
# reannotated_df.loc[reannotated_df.anno_lvl_0_final.isin(progenitors_labels),'anno_lvl_0_final'] = "PROGENITORS"
# reannotated_df.loc[reannotated_df.anno_lvl_0_final.isin(['DN(P) T CELL', 'DN(early) T CELL']),'anno_lvl_0_final'] = "NK/T CELLS"

# anno_groups_dict = {}
# for s in reannotated_df["anno_lvl_0_final"].unique():
#     anno_groups_dict[s] = reannotated_df[reannotated_df["anno_lvl_0_final"] == s]["anno_lvl_2_final_clean"].unique()


# In[59]:


def _plot_old_annotation_confusion(anno_group="B CELLS",
                                   min_cells = 5):
    keep_col = anno_groups_dict[anno_group]
    keep_row = annotations_conf_mat.index[np.any(annotations_conf_mat[anno_groups_dict[anno_group]] > min_cells, 1)]
    conf_mat = annotations_conf_mat.loc[keep_row, keep_col].T
    conf_mat = (conf_mat.T/conf_mat.sum(1)).T

    # row_order = conf_mat.max(1).sort_values(ascending=False).index
    row_order = conf_mat.index.sort_values()
    conf_mat = conf_mat.loc[row_order] ## Sort to have some sort of diagonal

    col_order = conf_mat.mean(0).sort_values(ascending=False).index
    conf_mat = conf_mat[col_order] ## Sort to have some sort of diagonal
    sns.heatmap(conf_mat,xticklabels=True, yticklabels=True)


# In[ ]:


plt.rcParams["figure.figsize"] = [20,12]
for g in anno_groups_dict.keys():
    _plot_old_annotation_confusion(g)
    plt.xlabel("old annotation");
    plt.ylabel("new annotation");
    plt.savefig(figdir + "oldVSnew_{g}.png".format(g="_".join(g.split("/"))), bbox_inches="tight")
    plt.show()


# In[24]:


anno_obs.to_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.clean.csv".format(t=timestamp))


# In[155]:


data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.annotated.clean.csv".format(t=timestamp)


# In[156]:


## Save file for lvl1 groupings 
anno_lvl_0_df = pd.DataFrame([(x,y) for x,y in anno_groups_dict_rev.items()])
anno_lvl_0_df.columns = ['anno_lvl_2_final_clean', "anno_lvl_0"]
anno_lvl_0_df.to_csv('/home/jovyan/mount/gdrive/Pan_fetal/annotations/annotation_terms_20211106.csv')


# ---
