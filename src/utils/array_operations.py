#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
"""
Created on Fri Nov 20 21:30:44 2020

@author: tengyao
"""

import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from collections.abc import Iterable

#%%
"""
match function finds positions of first match for elements in small_list in
    big_list.
small_list - a list, or numpy.ndarray, or pandas.Series object
big_list - a list, or numpy.ndarray, or pandas.Series object
nomatch - value to include if no match is found, -1 by default
[return] - a list of indices of first matches of small_list in big_list
"""
def match(small_list, big_list, nomatch=-1, sortable=True):
    if sortable:
        order = np.array(np.argsort(big_list))
        big_sorted = np.array(big_list)[order]
        small_list = np.array(small_list)
        l = np.searchsorted(big_sorted, small_list, side='left')
        insert_at_last = l == len(order)
        l[insert_at_last] = 0
        ifnomatch = insert_at_last | (big_sorted[l]!=small_list)
        ret = order[l]
        if np.any(ifnomatch):
            ret[ifnomatch] = nomatch
    else:
        ret = np.array([big_list.index(item) for item in small_list])
    return ret

#%%
"""
lookup function matches a vector of `lookup_value` to the `match_col` column of
    `dataframe` and returns the corresponding values in `result_col` colume of
    `dataframe`. 
lookup_value - a vector of values to lookup, can be a panda.Series, 
    numpy.ndarray or list
dataframe - the lookup table/array, can be a pandas.DataFrame, panda.Series or
    numpy.ndarray or list
match_col - column of `dataframe` to use as values to look up from, can be an 
    int or a string, if equal to -1, rownames of `dataframe` are used
result_col - column of `dataframe` to use as the result vector, can be an int
    or a string, if equal to -1, rownames of `dataframe` are used, if equals
    None, the matched ordering is returned.
[return] - a numpy.ndarray of lookup results, with unmatched positions filled
    with NaN.
"""
def lookup(lookup_value, dataframe, match_col=0, result_col=None):
    isIterable = isinstance(lookup_value, Iterable)
    lookup_value = pd.Series(lookup_value) if isIterable else pd.Series([lookup_value])
    dataframe = pd.DataFrame(dataframe)
    if type(match_col) is int:
        tmp = dataframe.iloc[:, match_col] if match_col >= 0 else dataframe.index
    else:
        tmp = dataframe[match_col]
       
    if result_col is None:
        ret = np.array(match(lookup_value, tmp))
        return(ret if isIterable else ret[0])
    elif type(result_col) is int:
        tmp2 = dataframe.iloc[:, result_col] if result_col >= 0 else dataframe.index
    else:
        tmp2 = dataframe[result_col]
    tmp2 = np.append(tmp2, np.nan)
    m = match(lookup_value, tmp, dataframe.shape[0])
    return(tmp2[m] if isIterable else tmp2[m][0])
    
#%%
"""
Slicing by exclusion
    vector: a list, numpy.array, pandas.Series or tuple object
    idx: indices to be excluded, can be a single integer or an Iterable of int
"""
def exclude(vector, idx):
    if not isinstance(idx, Iterable): idx = range(idx, idx+1)
    ret = [val for i, val in enumerate(vector) if i not in idx]
    if type(vector) is np.ndarray: return np.array(ret)
    if type(vector) is pd.core.series.Series: return pd.Series(ret)
    if type(vector) is tuple: return tuple(ret)
    return ret
    
#%%
"""
aggregate a matrix by a clustering of the column/rows
    X - a matrix in the form of an numpy.ndarray
    axis - axis to aggregate, 0 means aggregate by column, 1 mean aggregate by row
    func - function to aggregation, default is numpy.mean
return: a tuple of aggregated matrix X_agg (ndarray of dimension len(labels) in
    the axis direction) and a vector of labels
"""
X = np.array([[1,2,3,4],[4,5,6,8],[7,8,9,0]]).T
def aggregate(X, axis=0, by=None, func=np.mean):
    if by is None: by = np.zeros(X.shape[axis], int)
    by = np.array(by)
    labels = np.unique(by) # unique labels in the `by' vector
    if len(X.shape)==1:
        X_agg = np.ndarray(len(labels))
        for i, label in enumerate(labels):
            X_agg[i] = func(X[by==label])
    elif axis==1:
        X_agg = np.ndarray((X.shape[0], len(labels)))
        for i, label in enumerate(labels):
            X_agg[:, i] = func(X[:, by==label], axis=axis)
    else:
        X_agg = np.ndarray((len(labels), X.shape[1]))
        for i, label in enumerate(labels):
            X_agg[i, :] = func(X[by==label, :], axis=axis)
    return (X_agg, labels)
    
#%%
"""
select marker genes for each celltype
"""
def select_marker(X, by, genes, rankby='logFC', pval_cutoff=0.05, logFC_cutoff=0):    
    by = np.array(by)
    X_agg, unique_celltypes = aggregate(X, by=by)
    logFC = {}; max_pval = {}; best_celltype = {}
    for i, gene in enumerate(genes):
        order = np.argsort(X_agg[:, i])
        logFC[gene] = X_agg[order[-1], i] - X_agg[order[-2], i]
        best_celltype[gene] = unique_celltypes[order[-1]]
        tmp = 0
        for ct in exclude(unique_celltypes, order[-1]):
            _, pval = stats.ranksums(X[by==ct, i], X[by==best_celltype[gene], i])
            tmp = max(tmp, pval)
        max_pval[gene] = tmp
    df = pd.DataFrame({'logFC': logFC, 'max_pval': max_pval, 
                       'best_celltype': best_celltype})  # collect all genes and 
    df = df.iloc[np.logical_and(max_pval <= pval_cutoff,  logFC >= logFC_cutoff), :]
    df = df.iloc[np.array(np.argsort(df['max_pval'])), :]
    print(df)
    ret = {}
    for ct in unique_celltypes:
        print(ct)
        ret[ct] = df.iloc[np.array(df['best_celltype'])==ct, [0,1]]
    return(ret)
        
        
#%%
"""
Mann-Whitney-Wilcoxon test in a sparse matrix
"""
def wilcox_test_csc(mx, group1, show_progress=False):
    n, p = mx.shape # mx is a n x p matrix
    p_vals = np.zeros(p)
    logFC = np.zeros(p)
    np.seterr(divide='ignore') # ignore divide by zero warning in logFC calculation
    
    index_range = range(p)
    if show_progress:
        from tqdm import tqdm
        index_range = tqdm(index_range)
    for j in index_range:
        idx = mx[:, j].nonzero()  # location of nonzero entries in jth column
        val = np.array(mx[:, j][idx]).flatten() # values of nonzero entries
        idx = idx[0] # only keep row index values
        if len(idx) == 0:
            p_vals[j] = 1
            logFC[j] = 0
        else:
            num_zero = n - len(idx) # number of zeros in the jth column
            val = np.append(val, 0) # append 0 to values

            rank = sp.stats.rankdata(val) # compute rank of all values, including 0
            rank[rank > rank[-1]] += num_zero - 1 # adjust ranks of entries larger than 0 by ties in 0
            rank[-1] += (num_zero - 1) / 2 # adjust rank of 0 by ties

            # compute the rank sum of group1 elements
            idx_in_grp1 = np.in1d(idx, group1) # boolean vector of whether nonzero entry belongs to group1
            num_grp1 = len(group1) 
            num_nonzero_grp1 = np.sum(idx_in_grp1)
            num_zero_in_grp1 = num_grp1 - num_nonzero_grp1
            rank_sum = np.sum(rank[:-1][idx_in_grp1]) + rank[-1] * num_zero_in_grp1

            # compute the Z statistic and p-value
            mean = num_grp1 * (n+1) / 2
            sd = np.sqrt(num_grp1 * (n - num_grp1) / 12 * 
                         (n + 1 - (num_zero_in_grp1^3 - num_zero_in_grp1) / n / (n - 1))) # adjust standard deviation for ties
            Z_stat = (rank_sum - mean) / sd
            p_vals[j] = 1 - sp.stats.norm.cdf(Z_stat)
            
            # compute the log fold change
            sum1 = np.sum(val[:-1][idx_in_grp1])
            sum2 = np.sum(val) - sum1
            logFC[j] = np.log((sum1 / num_grp1) / (sum2 / (n - num_grp1)))
    return p_vals, logFC