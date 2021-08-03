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
import scipy

import sklearn

def _map_query_to_panfetal(
    query_adata,
    split,
    batch_name,
    scvi_outs_dir = "/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/",
    query_gene_id_col = "gene_ids",
    model_dir = None,
    timestamp = "20210429"
    ):
    '''
    Use scArches method to do online update of scVI model with query dataset
    
    Params:
    -------
    - query_adata: anndata object of query dataset to map
    - split: name of split to use (to chose the model)
    - batch_name: how do you want to label the new data you are adding? Used only if `query_adata.obs["bbk"]` is `None`
    - scvi_outs_dir: directory storing the scVI training outputs
    - query_gene_id_col: string indicating column in `query_adata.var` containing geneIDs 
    (used for feature matching if the saved features from trained model are not matched)
    - model_dir: which model directory in `scvi_outs` to use? (Default: None, uses `scvi_{split}_model/`)
    
    Outputs:
    - anndata with new embedding in obsm ("X_scVI_project")
    - vae: model
    '''
    if model_dir is None:
        model_dir='scvi_' + split + '_model/'
    ## Check that the feature used for reference training match the var_names in query data
    var_names_model = pd.read_csv(scvi_outs_dir + model_dir + "var_names.csv", header=None)[0].values

    vars_overlap = any(pd.Series(var_names_model).isin(query_adata.var_names))

    ## Extract gene ids if not
    if not vars_overlap:
        raise ValueError("Vars don't overlap -- convert to geneIDs")
#         print("Vars don't overlap -- converting to geneIDs")
#         adata_ref_var = pd.read_csv(scvi_outs_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.var.csv'.format(t=timestamp, s=split), index_col=0)
#         adata_ref_var.iloc[var_names_model]['GeneID'].to_csv(scvi_outs_dir + model_dir + "var_names.csv", header=None, index=False)
#         query_adata.var["gene_name"] = query_adata.var_names.copy()
#         query_adata.var_names = query_adata.var[query_gene_id_col]
#         var_names_model = pd.read_csv(scvi_outs_dir + model_dir + "var_names.csv", header=None)[0].values
    
    ## Add batch column
    if "bbk" not in query_adata.obs.columns:
        print("considering all query cells from the same technical batch")
        query_adata.obs["bbk"] = batch_name
        
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
    vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
    query_adata_filled.obsm["X_scvi"] = vae_q.get_latent_representation()
    return(query_adata_filled, vae_q)

## I/O utils ##
def merge_query_and_reference(
    query_h5ad_file,
    split,
    ref_data_dir = '/home/jupyter/mount/gdrive/Pan_fetal/data4gpu_node/',
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

    concat_adata = ref_adata.concatenate(query_adata_mapped, batch_key="dataset", batch_categories=["reference", "query"], join="outer")
    sc.pp.neighbors(concat_adata, use_rep="X_scvi", n_neighbors = 30)
    return(concat_adata)

### Utils to propagate cell type labels ### 

def predict_label(adata, anno_col, 
                  neighbors_key = 'connectivities',
                  min_score = 0
                 ):
    min_ref_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]/10
    missing_anno = adata.obs["dataset"] == "query"

    ## Find neighbors of cells with conflicting annotation
    knn_graph = adata.obsp[neighbors_key]
    knn_graph_query = knn_graph[missing_anno,:]
    knn_graph_query[knn_graph_query.nonzero()] = 1

    ## Find most abundant cell label in neighbors
    neighbors_labels = pd.DataFrame()
    n_neighbors_labels = pd.DataFrame()

    annos = adata.obs[anno_col].copy()

    dummy_df = pd.get_dummies(annos)
    dummy_mat = scipy.sparse.csr_matrix(dummy_df)

    new_anno = knn_graph_query.dot(dummy_mat).toarray()

    n_neighbors = np.array(knn_graph_query.sum(1)).flatten()
    n_neighbors_ref = new_anno.sum(axis=1)
    new_anno_prob = new_anno.T/n_neighbors_ref
    new_anno_prob[np.isnan(new_anno_prob)] = 0
    new_anno_prob[:,n_neighbors_ref < min_ref_neighbors] = 0

    best_label = dummy_df.columns[new_anno_prob.argmax(0)].values
    best_label_score = new_anno_prob.max(0)
    
    best_label[best_label_score <= min_score] = "low_confidence"
    
    adata.obs.loc[missing_anno,'predicted_anno'] = best_label
    adata.obs.loc[missing_anno,'predicted_anno_prob'] = best_label_score

def plot_confusion_mat(adata, query_anno_col, show_low_confidence=True, **kwargs):
    missing_anno = adata.obs["dataset"] == "query"
    if show_low_confidence:
        conf_mat = sc.metrics.confusion_matrix(query_anno_col,"predicted_anno",  adata[missing_anno].obs, normalize=True)
    else:
        high_conf_prediction = adata.obs["predicted_anno"] != "low_confidence"
        conf_mat = sc.metrics.confusion_matrix(query_anno_col,"predicted_anno",  adata[missing_anno & high_conf_prediction].obs, normalize=True)
#     conf_mat.columns = conf_mat.columns.str.upper()
#     conf_mat.columns = conf_mat.columns.str.replace(" ", "_")
    col_order = conf_mat.idxmax(0).sort_values().index
    conf_mat = conf_mat[col_order] ## Sort to have some sort of diagonal
    sns.heatmap(conf_mat, xticklabels=True, yticklabels=True, **kwargs);
    plt.xlabel("Predicted label");
    plt.ylabel("Query label");

def plot_predicted_anno_probability(adata):
    missing_anno = adata.obs["dataset"] == "query"
    best_label_df = adata.obs[missing_anno]

    sns.stripplot(data=best_label_df, y="predicted_anno", x="predicted_anno_prob", orient='h', color="black");
    sns.boxplot(data=best_label_df, y="predicted_anno", x="predicted_anno_prob", orient='h', color=None);
    plt.xlabel("probability")
    
### Evaluation metrics ### 
def compute_silhouette(query_adata_mapped, query_anno_col = "annotation_V2"):
    sil_score_cells = sklearn.metrics.silhouette_samples(query_adata_mapped.obsm["X_scvi"], query_adata_mapped.obs[query_anno_col])
    query_adata_mapped.obs["sil_score"] = sil_score_cells
    sil_df = query_adata_mapped.obs[[query_anno_col, "sil_score"]].groupby(query_anno_col).mean()
    return(sil_df)

def _anno_silhouette_permutation(merged_adata, query_anno_col, n_permutations=1000):
    query_adata_mapped = merged_adata[merged_adata.obs["dataset"] == "query"]
    query_adata_mapped.obs[query_anno_col + "_random"] = query_adata_mapped.obs[query_anno_col].copy()
    rand_sil_df = pd.DataFrame()
    for i in range(n_permutations):
        np.random.shuffle(query_adata_mapped.obs[query_anno_col + "_random"].values)
        rand_sil_df = pd.concat([rand_sil_df, compute_silhouette(query_adata_mapped, query_anno_col + "_random")],1 )

    real_sil_df = compute_silhouette(query_adata_mapped, query_anno_col)
    for i in range(real_sil_df.shape[0]):
        real_sil_df["pval"] = 0
        real_sil_df.iloc[i,1] = sum(rand_sil_df.iloc[0] >= real_sil_df.iloc[0][0])/n_permutations
    return(real_sil_df)

## Calculate normalized mutual information
def compute_nmi(query_adata_mapped, query_anno_col):
    keep_ixs = query_adata_mapped.obs["predicted_anno"] != "low_confidence"
    nmi=sklearn.metrics.normalized_mutual_info_score(query_adata_mapped.obs[keep_ixs][query_anno_col], query_adata_mapped.obs[keep_ixs]["predicted_anno"])
    return(nmi)