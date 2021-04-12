import os,sys
import numpy as np 
import pandas as pd
import scanpy as sc
from packaging import version

import pandas as pd
import numpy as np
import scipy
import sys
from annoy import AnnoyIndex
from packaging import version
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from umap.umap_ import fuzzy_simplicial_set
from sklearn.neighbors import KDTree
from sklearn.neighbors import DistanceMetric
from sklearn.linear_model import Ridge
try:
	from scanpy import logging as logg
except ImportError:
	pass
try:
	import anndata
except ImportError:
	pass
try:
	import faiss
except ImportError:
	pass
from bbknn import * 

def get_graph_directed(pca, batch_list, ordered_batches, sorting_codes,tot_neighbors, n_pcs, approx, metric, use_faiss, n_trees):
    #in case we're gonna be faissing, turn the data to float32
    if metric=='euclidean' and not approx and 'faiss' in sys.modules and use_faiss:
        pca = pca.astype('float32')
    # knn_indices_shape = max([sum(max_neighbors_within_batch - np.abs(sorting_codes - sorting_codes[i])) for i in range(len(batches))])    
    ## Make array of n_neighbors for each batch pair
    knn_batches = np.zeros((len(ordered_batches), len(ordered_batches))).astype(int)
    unique, counts = np.unique(batch_list, return_counts=True)
    ## Assign k to each batch based on distance in time (must sum to tot_neighbors)
    for from_ind in range(len(ordered_batches)):
        time_sim = 1 / (1 + np.abs(sorting_codes - sorting_codes[from_ind]))
        k_prob = time_sim/sum(time_sim)
        batch_neighbors = np.round(k_prob*tot_neighbors)
        ## Add missing neighbors to highest ks to get to max
        error = int(tot_neighbors - sum(batch_neighbors))
        if error > 0:
            batch_neighbors[batch_neighbors.argsort()[-error:]] = batch_neighbors[batch_neighbors.argsort()[-error:]] + 1 
        if error < 0:
            batch_neighbors[batch_neighbors.argsort()[error:]] = batch_neighbors[batch_neighbors.argsort()[error:]] - 1 
        ## check that all the ordered_batches have sufficient number of cells
        if any(batch_neighbors > counts):
            raise ValueError("Not all ordered_batches have at least `batch_neighbors` cells in them.")
        knn_batches[from_ind,:] = batch_neighbors
    #create the output matrices, with the indices as integers and distances as floats
    knn_distances = np.zeros((pca.shape[0], tot_neighbors))
    knn_indices = np.copy(knn_distances).astype(int)
    # Initialize column range to 0 for each batch  
    col_range_start = np.zeros(len(ordered_batches)).astype("int")
    for to_ind in range(len(ordered_batches)):
        #this is the batch that will be used as the neighbour pool
        #create a boolean mask identifying the cells within this batch
        #and then get the corresponding row numbers for later use
        batch_to = ordered_batches[to_ind]
        batch_list = np.asarray([str(i) for i in batch_list])
        mask_to = batch_list == batch_to
        ind_to = np.arange(len(batch_list))[mask_to]
         #create the faiss/cKDTree/KDTree/annoy, depending on approx/metric
        ckd = create_tree(data=pca[mask_to,:n_pcs],approx=approx,metric=metric,
                          use_faiss=use_faiss,n_trees=n_trees)
        for from_ind in range(len(ordered_batches)):   
            #this is the batch that will have its neighbours identified
            #repeat the mask/row number getting
            batch_from = ordered_batches[from_ind]
            mask_from = batch_list == batch_from
            ind_from = np.arange(len(batch_list))[mask_from]
            neighbors_within_batch = knn_batches[from_ind,to_ind]
            if neighbors_within_batch > 0:
                ckdout = query_tree(data=pca[mask_from,:n_pcs],ckd=ckd,
                                    neighbors_within_batch=neighbors_within_batch,
                                    approx=approx,metric=metric,use_faiss=use_faiss)
                #the identified indices are relative to the subsetted PCA matrix
                #so we need to convert it back to the original row numbers
                for i in range(ckdout[1].shape[0]):
                    for j in range(ckdout[1].shape[1]):
                        ckdout[1][i,j] = ind_to[ckdout[1][i,j]]
                #save the results within the appropriate rows and columns of the structures
                col_range = np.arange(col_range_start[from_ind], col_range_start[from_ind] + neighbors_within_batch)
                knn_indices[ind_from[:,None],col_range[None,:]] = ckdout[1]
                knn_distances[ind_from[:,None],col_range[None,:]] = ckdout[0]
                col_range_start[from_ind] = col_range_start[from_ind] + neighbors_within_batch
    return(knn_distances, knn_indices)

def directed_bbknn_pca_matrix(pca, batch_list, ordered_batches, sorting_codes,tot_neighbors=30, 
                              n_pcs=30, trim=None,
                              approx=True, n_trees=10, use_faiss=True, metric='angular',
                              set_op_mix_ratio=1, local_connectivity=1):
    '''
	BBKNN variant that runs on a PCA matrix and list of per-cell batch assignments instead of
	an AnnData object. While regular BBKNN searches for an equal number of nearest neighbors between all batches, 
    directed BBKNN takes into account a user defined ordering of batches, searching more nearest neighbors between
    batches that are closer in order (e.g. closer time-point). 
	Returns a ``(distances, connectivities, parameters)`` tuple, like what would have been stored in the AnnData object.
	The connectivities are the actual neighbourhood graph.
	Input
	-----
	pca : ``numpy.array``
		PCA (or other dimensionality reduction) coordinates for each cell, with cells as rows.
	batch_list : ``numpy.array`` or ``list``
		A list of batch assignments for each cell.
	'''
    if pca.shape[0] != len(batch_list):
        raise ValueError("Different cell counts indicated by `pca.shape[0]` and `len(batch_list)`.")
    #convert batch_list to np.array of strings for ease of mask making later
    batch_list = np.asarray([str(i) for i in batch_list])
    knn_distances, knn_indices = get_graph_directed(pca=pca,batch_list=batch_list,ordered_batches=ordered_batches,sorting_codes=sorting_codes,tot_neighbors=tot_neighbors,
                                                    n_pcs=n_pcs,n_trees=n_trees,approx=approx,metric=metric,use_faiss=use_faiss)
        #sort the neighbours so that they're actually in order from closest to furthest
    newidx = np.argsort(knn_distances,axis=1)
    knn_indices = knn_indices[np.arange(np.shape(knn_indices)[0])[:,np.newaxis],newidx]
    knn_distances = knn_distances[np.arange(np.shape(knn_distances)[0])[:,np.newaxis],newidx]

    #this part of the processing is akin to scanpy.api.neighbors()
    dist, cnts = compute_connectivities_umap(knn_indices, knn_distances, knn_indices.shape[0],
                                             knn_indices.shape[1], set_op_mix_ratio=set_op_mix_ratio,
                                             local_connectivity=local_connectivity)

    #trimming. compute default range if absent
    if trim is None:
        trim = 10 * knn_distances.shape[1]
    #skip trimming if set to 0, otherwise trim
    if trim > 0:
        cnts = trimming(cnts=cnts,trim=trim)
    #create a collated parameters dictionary
    #determine which neighbour computation was used, mirroring create_tree() logic
    if approx:
        computation='annoy'
    elif metric == 'euclidean':
        if 'faiss' in sys.modules and use_faiss:
            computation='faiss'
        else:
            computation='cKDTree'
    else:
        computation='KDTree'
    #we'll have a zero distance for our cell of origin, and nonzero for every other neighbour computed
    params = {'n_neighbors': len(dist[0,:].data)+1, 'method': 'umap', 
              'metric': metric, 'n_pcs': n_pcs, 
              'bbknn': {'trim': trim, 'computation': computation}}
    return(dist, cnts, params)

def directed_bbknn(adata, batch_key='batch', sorting_key = "age", tot_neighbors=30, use_rep='X_pca', approx=True, metric='angular', copy=False, **kwargs):
    sorted_df = adata.obs[[batch_key, sorting_key]] \
        .drop_duplicates() \
        .sort_values(sorting_key, na_position="first") 
    ## Assign order index to batch list 
    sorting_codes = sorted_df[sorting_key].astype("category").values.codes + 1
    ordered_batches = sorted_df[batch_key].values
    batch_list = adata.obs[batch_key]
    #prepare bbknn_pca_matrix input
    pca = adata.obsm[use_rep]
    ## Call directed BBKNN
    bbknn_out = directed_bbknn_pca_matrix(pca=pca, batch_list=batch_list,
                                          ordered_batches=ordered_batches,
                                          sorting_codes=sorting_codes,
                                          tot_neighbors=tot_neighbors,
                                 approx=approx, metric=metric, **kwargs)
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = bbknn_out[2]
    adata.uns['neighbors']['params']['use_rep'] = use_rep
    adata.uns['neighbors']['params']['bbknn']['batch_key'] = batch_key
    #store the graphs in .uns['neighbors'] or .obsp, conditional on anndata version
#     if version.parse(str(anndata.__version__)) < version.parse('0.7.0'):
#         adata.uns['neighbors']['distances'] = bbknn_out[0]
#         adata.uns['neighbors']['connectivities'] = bbknn_out[1]
#         logg.info('	finished', time=start,
#             deep=('added to `.uns[\'neighbors\']`\n'
#             '	\'distances\', distances for each pair of neighbors\n'
#             '	\'connectivities\', weighted adjacency matrix'))
#     else:
    adata.obsp['distances'] = bbknn_out[0]
    adata.obsp['connectivities'] = bbknn_out[1]
    adata.uns['neighbors']['distances_key'] = 'distances'
    adata.uns['neighbors']['connectivities_key'] = 'connectivities'
    return adata if copy else None