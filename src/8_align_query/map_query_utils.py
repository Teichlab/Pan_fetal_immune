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


# def _add_all_query_genes(merged_adata, query_adata_full):
#     if not any(merged_adata.var_names.isin(query_adata_full.var_names)):
#         raise ValueError("var_names don't match between query and merged AnnData")

#     if not any(query_adata_full.obs_names.str.endswith("-query")):
#         query_adata_full.obs_names = query_adata_full.obs_names + "-query"

#     ## Do the merge
#     full_merged_adata = anndata.concat([merged_adata, query_adata_full[0:142]], axis=1, join="outer", merge="unique", uns_merge="unique")

#     ## Check that the number of obs is right
#     if not full_merged_adata.n_obs == merged_adata.n_obs:
#         raise AssertionError("The number of obs doesn't match, something is wrong in your join")

#     ## Check that you have more expression vals than before
#     one_query_cell = query_adata_full.obs_names[10]
#     if not len(full_merged_adata[one_query_cell].X.nonzero()[0]) > len(merged_adata[one_query_cell].X.nonzero()[0]):
#         raise AssertionError("You have less or the same expression values for query cells than before. Are you sure that query_adata_full is the dataset before feature selection?")

#     return(full_merged_adata)

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

## Functions to compute MNN/KNN similarity ratio 
def _find_MNNs(merged_adata, k=30, n_jobs=5):
    '''
    Find mutual nearest neighbors between query and reference data
    '''
    from scipy.spatial import cKDTree

    ## Extract embedding
    X_emb = merged_adata.obsm["X_scvi"].copy()

    is_query = merged_adata.obs["dataset"] == "query"
    is_reference = merged_adata.obs["dataset"] == "reference"

    X_emb_ref = X_emb[is_reference,:]
    X_emb_que = X_emb[is_query,:]

    ## Find mutual nearest neighbors pairs
    # ## Borrowed from chriscainx/mnnpy
    k1=k2=k
    data_query = X_emb_que
    data_ref = X_emb_ref
    k_index_ref = cKDTree(data_ref).query(x=data_query, k=k1, n_jobs=n_jobs)[1]

    ## Subset to reference cells that have some nn
    ref_w_nn = np.unique(k_index_ref.flatten())
    k_index_query = cKDTree(data_query).query(x=data_ref[ref_w_nn,:], k=k2, n_jobs=n_jobs)[1]
    mutual_ref = []
    mutual_query = []
    for index_ref in range(len(ref_w_nn)):
        for index_query in k_index_query[index_ref]:
            if ref_w_nn[index_ref] in k_index_ref[index_query]:
                mutual_ref.append(ref_w_nn[index_ref])
                mutual_query.append(index_query)

    mutual_query = np.array(mutual_query)
    mutual_ref = np.array(mutual_ref)

    ## Convert to array of MNNs
    ref_obs = merged_adata.obs_names[is_reference]
    que_obs = merged_adata.obs_names[is_query]
    #     mnn_reference = np.zeros(shape=(ref_obs.shape[0], k))
    mnn_reference = np.empty(shape=(ref_obs.shape[0], k))
    mnn_reference[:] = np.nan

    has_mnn_ixs = np.unique(mutual_ref)

    for i in has_mnn_ixs:
        ds = mutual_query[mutual_ref==i]
        mnn_reference[i,0:ds.shape[0]] = ds.copy()

    mnn_query = np.empty(shape=(que_obs.shape[0], k))
    mnn_query[:] = np.nan
    has_mnn_ixs = np.unique(mutual_query)
    for i in has_mnn_ixs:
        ds = mutual_ref[mutual_query==i]
        mnn_query[i,0:ds.shape[0]] = ds.copy()

    mnn_reference[~np.isnan(mnn_reference)] = mnn_reference[~np.isnan(mnn_reference)].astype("int")
    mnn_query[~np.isnan(mnn_query)] = mnn_query[~np.isnan(mnn_query)].astype("int")

    return(mnn_query, mnn_reference)

def _scArches_adjusted_dist(dist_nns):
    # compute standard deviation
    std_d = np.sqrt(sum(np.sqrt(dist_nns))/dist_nns.shape[0])
    # apply gaussian kernel
    adj_dist_nns = np.exp(-dist_nns/np.square(2/std_d))
    return(adj_dist_nns)

def _MNN_to_KNN_similarity_ratio(merged_adata, mnn_reference, mnn_query):
    from pynndescent import NNDescent
    ## Extract embedding
    X_emb = merged_adata.obsm["X_scvi"].copy()

    is_query = merged_adata.obs["dataset"] == "query"
    is_reference = merged_adata.obs["dataset"] == "reference"

    has_mnn = np.unique(mnn_query[~np.isnan(mnn_query)].astype("int"))

    ## Mean similarity between mutual nearest neighbors
    dists = scipy.spatial.distance.cdist(X_emb[is_query,:], X_emb[is_reference,:][has_mnn,:], metric="euclidean")

    mnn_dists_reference = mnn_reference.astype("float64").copy()
    for i in range(len(has_mnn)):
        mnns_ixs = mnn_reference[has_mnn[i]][~np.isnan(mnn_reference[has_mnn[i],:])].astype("int").ravel()
        mnn_dists_reference[has_mnn[i],~np.isnan(mnn_reference[has_mnn[i],:])] = dists[mnns_ixs,i]

    mnn_dists_query = mnn_query.astype("float64").copy()
    for i in range(sum(is_query)):
        mnns_ixs = mnn_query[i, :][~np.isnan(mnn_query[i,:])].ravel()
        for j in range(len(mnns_ixs)): 
            mnn_dists_query[i,np.where(~np.isnan(mnn_query[i,:]))[0][j]] = dists[i,np.where(has_mnn==mnns_ixs[j])[0]]

    ## Adjust distances as in scArches paper
    adj_mnn_dists_reference = mnn_dists_reference.copy()
    for i in range(len(mnn_dists_reference)):
        if adj_mnn_dists_reference[i,mnn_dists_reference[i] > 0].shape[0] > 0:
            adj_mnn_dists_reference[i,mnn_dists_reference[i] > 0] = _scArches_adjusted_dist(mnn_dists_reference[i, mnn_dists_reference[i] > 0])

    adj_mnn_dists_query = mnn_dists_query.copy()
    for i in range(len(mnn_dists_query)):
        if adj_mnn_dists_query[i,mnn_dists_query[i] > 0].shape[0] > 0:
            adj_mnn_dists_query[i,mnn_dists_query[i] > 0] = _scArches_adjusted_dist(mnn_dists_query[i, mnn_dists_query[i] > 0])

    adj_mnn_dists_query[np.isnan(adj_mnn_dists_query)] = 0
    adj_mnn_dists_reference[np.isnan(adj_mnn_dists_reference)] = 0

    mean_mnn_sim_reference = adj_mnn_dists_reference.mean(1)
    mean_mnn_sim_query = adj_mnn_dists_query.mean(1)

    ## Compute distances between nearest neighbors
    knns = NNDescent(X_emb, metric="euclidean").query(X_emb, k=30)
    knn_dists = knns[1][:,1:]

    adj_knn_dists = knn_dists.copy()
    for i in range(len(knn_dists)):
        adj_knn_dists[i] = _scArches_adjusted_dist(knn_dists[i])

    ## Mean similarity between nearest neighbors
    mean_knn_sim_reference = adj_knn_dists[is_reference].mean(1)
    mean_knn_sim_query = adj_knn_dists[is_query].mean(1)

    mnn_sim_ratio_reference = mean_mnn_sim_reference/mean_knn_sim_reference
    mnn_sim_ratio_query = mean_mnn_sim_query/mean_knn_sim_query

    ## Save everything in adata
    merged_adata.obs['knn_sim'] = 0
    merged_adata.obs.loc[is_reference,'knn_sim'] = mean_knn_sim_reference
    merged_adata.obs.loc[is_query,'knn_sim'] = mean_knn_sim_query

    merged_adata.obs['mnn_sim'] = 0
    merged_adata.obs.loc[is_reference,'mnn_sim'] = mean_mnn_sim_reference
    merged_adata.obs.loc[is_query,'mnn_sim'] = mean_mnn_sim_query

    merged_adata.obs['mnn_sim_ratio'] = 0
    merged_adata.obs.loc[is_reference,'mnn_sim_ratio'] = mnn_sim_ratio_reference
    merged_adata.obs.loc[is_query,'mnn_sim_ratio'] = mnn_sim_ratio_query