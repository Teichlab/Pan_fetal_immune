### Utils ###
import sys,os
import pandas as pd
import scanpy as sc
import numpy as np
import scipy
from sklearn.linear_model import Ridge
import time
from bbknn import bbknn

sys.path.append('.')
from genes import cc_genes
from genes import IG_genes

### Functions used for preprocessing ###

def remove_geneset(adata,geneset):
    ## Check that var_names are not ensemblIDs
    # if adata.var_names.str.startswith('ENS')[0]:
    #     raise ValueError('adata.var_names are ensemblIDs')
    adata = adata[:,~adata.var_names.isin(list(geneset))].copy()
    return adata

def is_cycling(adata,cc_genes=cc_genes,cut_off=0.4):
    X = np.mean(adata.X[:,adata.var_names.isin(cc_genes)],axis=1)
    plt.hist(X)
    adata.obs['Cycle_score'] = X
    adata.obs['isCycle'] = X>cut_off
    
def pfi_preprocess(sdata, how="pd", inplace=False, genesets2remove=[cc_genes]):
    '''
    General preprocessing function
    
    how options:
    - 'p': run pca
    - 'd': run diffusion map
    '''
    if inplace:
        sdata_pp = sdata
    else:
        sdata_pp = sdata.copy()
    start = time.time()
    sc.pp.highly_variable_genes(sdata_pp, min_mean=0.001, max_mean=10, subset=True)
    sc.pp.scale(sdata_pp,max_value=10)
    for s in genesets2remove:
        sdata_pp = remove_geneset(sdata_pp,s)
    if "p" in how:
        sc.pp.pca(sdata_pp, use_highly_variable=True)
    if "d" in how:
        sc.pp.neighbors(sdata_pp)
        sc.tl.diffmap(sdata_pp)
    pp_time = time.time() - start
    print("Preprocessing runtime: ", str(pp_time))
    return(sdata_pp)

def pfi_clustering(adata, how="pbul", batch_key = "bbk",
                   res=0.5,
                   use_highly_variable=True, plot=True):
    if "p" in how:
        sc.pp.pca(adata, use_highly_variable=use_highly_variable)
        if plot:
            sc.pl.pca(adata, color="method", components=["1,2", '3,4', '5,6'])
    if "b" in how:
        start=time.time()
        bbknn(adata, batch_key = batch_key, n_pcs=30, approx=True)
        bbknn_time = time.time()-start
        print("BBKNN runtime: ", str(bbknn_time))
    if 'u' in how:
        start=time.time()
        sc.tl.umap(adata)
        umap_time = time.time()-start
        print("UMAP runtime: ", str(umap_time))
        if plot:
            sc.pl.umap(adata, color=["method", batch_key])
    if 'l' in how:
        lab_ls = str(res).split('.')
        if len(lab_ls)==2:
            label = ''.join(lab_ls) + "0"
        else:
            label = ''.join(lab_ls)
        start=time.time()
        sc.tl.leiden(adata, resolution=res, key_added='leiden_' + label, n_iterations=5)
        cl_time = time.time()-start
        print("Leiden clustering runtime: ", str(cl_time))
        
def ridge_regression(adata,batch_key,confounder_key=[], chunksize=1e8):
	'''
	batch regression tool (improved for efficiency by K. Polanski)
    ---
	batch_key = list of observation categories to be regressed out
	confounder_key = list of observation categories to be kept
	chunksize = how many elements of X to process at once, will iterate over genes
	---
    output: adds 'X_explained' and 'X_remain' layers to adata
	'''
	
	dummy = pd.get_dummies(adata.obs[batch_key+confounder_key],drop_first=False)
	if len(batch_key)>1:
		batch_index = np.logical_or.reduce(np.vstack([dummy.columns.str.startswith(x) for x in batch_key]))
	else:
		batch_index = np.vstack([dummy.columns.str.startswith(x) for x in batch_key])[0]
	dm = np.array(dummy)[:,batch_index]
	
	LR = Ridge(fit_intercept=False, alpha=1.0)
	chunkcount = np.ceil(chunksize/adata.shape[0])
	X_explained = []
	X_remain = []
	for ind in np.arange(0,adata.shape[1],chunkcount):
		X_exp = adata.X[:,np.int(ind):np.int(ind+chunkcount)] # scaled data
		if scipy.sparse.issparse(X_exp):
			X_exp = X_exp.todense()
		LR.fit(dummy,X_exp)	
		X_explained.append(dm.dot(LR.coef_[:,batch_index].T))
		X_remain.append(X_exp - X_explained[-1])
	
	X_explained = np.hstack(X_explained)
	X_remain = np.hstack(X_remain)
	adata.layers['X_remain'] = X_remain
	adata.layers['X_explained'] = X_explained
    
    
def _propagate_labels(adata, anno_col):
    '''
    Propagate labels to nan cells based on KNN graph
    '''
    anno_nans = (adata.obs[anno_col]=="nan").values
    nan2labelled_conns = adata.obsp["connectivities"][:,anno_nans]
    ## Get KNN edges between nans and all cells
    nan2labelled_conns[nan2labelled_conns.nonzero()] = 1
    
    ## Make dummy matrix of labels
    lab_dummies = pd.get_dummies(adata.obs[anno_col])
    lab_unique = lab_dummies.columns
    lab = lab_dummies.to_numpy().T
    
    ## Calculate label probability based on labels of neighbours
    class_prob =  lab @ nan2labelled_conns
    norm = np.linalg.norm(class_prob, 2, axis=0)
    class_prob = class_prob / norm
    class_prob = (class_prob.T - class_prob.min(1)) / class_prob.ptp(1)
    
    ## Pick label with max prob
    new_labs = lab_unique[class_prob[:,lab_unique!="nan"].argmax(1)].to_numpy()
    max_prob = class_prob[:,lab_unique!="nan"].max(1)
    new_labs[max_prob==0] = np.nan # exclude cells with no neighbour != nan

    adata.obs[anno_col + "_propagated"] = adata.obs[anno_col]
    adata.obs[anno_col + "_propagated"].loc[anno_nans] = new_labs
    
    
### I/O utils ###

def _load_split_and_annotation(split, full_mat=False, PFI_prefix = "PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118", data_dir = '/nfs/team205/ed6/data/Fetal_immune/'):
    '''
    Load anndata of Pan Fetal Immune split + annotation for the same split
    '''
    if not full_mat:
        split_file = "{prefix}.{split}.batchCorrected.h5ad".format(prefix = PFI_prefix, split = split)
    else:
        split_file = "{prefix}.{split}.h5ad".format(prefix = PFI_prefix, split = split)
    adata = sc.read_h5ad(data_dir + split_file)
    adata_obs = pd.read_csv("/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.wGut.full_obs.annotated.csv", index_col=0)

    ### Load manual annotations
    anno_dir = '/nfs/team205/ed6/bin/Pan_fetal_immune/manual_annotation/'
    keep_anno = [split]

    adata_obs["anno_lvl_1"] = np.nan
    adata_obs["anno_lvl_2"] = np.nan

    for anno in keep_anno:
        anno_file = "{prefix}.{anno}.batchCorrected_annotation.csv".format(prefix = PFI_prefix, anno = anno)
        if anno_file in os.listdir(anno_dir):
            anno_df = pd.read_csv(anno_dir + anno_file, index_col=0)
        else:
            print("No .csv found for annotation" + anno)
        ## Check for collisions between annotations (they should have been manually fixed)
        if adata_obs.loc[anno_df.index,'anno_lvl_1'].isna().all():        
            adata_obs.loc[anno_df.index,"anno_lvl_1"] = anno_df["anno_lvl_1"]
            adata_obs.loc[anno_df.index,"anno_lvl_2"] = anno_df["anno_lvl_2"]
        else:
            n_cells=anno_df.index.shape[0] - adata_obs.loc[anno_df.index,'anno_lvl_1'].isna().sum()
            print("Error! {n} cells are already annotated".format(n=n_cells))

    adata.obs = pd.concat([adata.obs, adata_obs.loc[adata.obs_names][["anno_lvl_1", "anno_lvl_2"]]], 1)
    return(adata)