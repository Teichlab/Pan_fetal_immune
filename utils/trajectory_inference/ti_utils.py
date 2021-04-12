import matplotlib
import seaborn as sns
import scipy.sparse
import scanpy as sc
import pandas as pd
import numpy as np
import cellrank as cr

def summarise_connectivities(adata, obs_col, plot=True, key="connectivities"):
    '''
    Make a matrix summarizing transitions between cells with the same value of obs column
    (Useful to visualize directed transition matrices)
    '''
    ## From connectivities to adjacency
    conn_biased_adj = adata.obsp[key]
    ## transform obs columns to dummy variables
    obs_dummies = pd.get_dummies(adata.obs[obs_col])
    obs_dummies_mat = scipy.sparse.csr_matrix(obs_dummies.values)
    ## average connections over obs
    conn_obs = obs_dummies_mat.T.dot(conn_biased_adj).dot(obs_dummies_mat).toarray()
    conn_obs = pd.DataFrame(conn_obs, columns=obs_dummies.columns, index=obs_dummies.columns)
    if plot:
        sns.heatmap(conn_obs)
    else:
        return(conn_obs)
   
### Cellrank analysis

def make_kernel(adata, kernel_type='ac', weights={'c':0.8, "a":0.2}, pk_k=3):
    '''
    Make kernel for cellrank analysis with combo of kernels
    
    - adata: the anndata
    - kernel_type: 'a' for PalantirKernel using "age" as time, 'c' for ConnectivityKernel
    - weights: should be the same length as type, weights to give to each kernel
    '''
    if "a" in kernel_type:
        adata.obs["age"] = adata.obs["age"].astype('int64')
        age_k = cr.tl.kernels.PalantirKernel(adata, time_key="age").compute_transition_matrix(k=pk_k)
        kern = age_k
    if "c" in kernel_type:
        conn_k = cr.tl.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        kern = conn_k
    if "c" in kernel_type and "a" in kernel_type:
        combo_k = weights["c"]*conn_k + weights["a"]*age_k
        kern = combo_k
    return(kern)
    
    
def run_cellrank(adata, kernel_type="ac", weights={'c':0.8, "a":0.2}, compute_macrostates=False):
    ## Convert to category to use to name macrostates
    adata.obs["anno_lvl_2"] = adata.obs["anno_lvl_2"].astype("category")

    ## Use the number of annotated clusters as number of macrostates/schur vectors
    n_states = adata.obs["anno_lvl_2"].unique().shape[0]
    
    ## Build kernel and estimator
    kern = make_kernel(adata, kernel_type = kernel_type, weights= weights)
    g = cr.tl.estimators.GPCCA(kern)
    
    ## Compute schur vectors
    g.compute_schur(n_components=n_states, method="krylov")
    
    ## Estimate generalized pseudotime
    
    if compute_macrostates:
        g.compute_macrostates(n_states=n_states, cluster_key="anno_lvl_2", n_cells=30)
        g.set_terminal_states_from_macrostates()
    