from collections import Counter
from collections import defaultdict
import scanpy as sc
import scrublet as scr
import pandas as pd
import pickle as pkl
from .markers import find_markers, show_marker_plot
from .model import generate_training_X
import numpy as np
from bbknn import bbknn

import scipy
import matplotlib.pyplot as plt
import seaborn as sns
import re
import glob
import os
import sys
from geosketch import gs

def create_directory(folder_name):
    if not os.path.exists(folder_name):
        print('no existing folder, making {}...'.format(folder_name))
        print(os.system('mkdir {}'.format(folder_name)))
    else:
        print('{} folder already exists'.format(folder_name))

def calculate_cell_count(adata,anno_key,donor_key = 'donor'):

    meta_table = adata.obs
    donor_list = np.unique(meta_table[donor_key])
    count_dict = {}

    # uncorrected raw counts
    s_list = []
    for donor in donor_list:
        c_donor = (meta_table[donor_key] == donor)
        if np.sum(c_donor)==0:
            continue
        order = str(int(adata.obs['Order'][c_donor][0]))
        age = adata.obs['Age'][c_donor][0]
        stage = adata.obs['Stage'][c_donor][0]

        s_list.append(meta_table[c_donor].groupby(anno_key)['Sample'].count())
        s_list[-1] = s_list[-1].rename("%s_%s_%s_%s"%(order, donor, age, stage))

    s_list = sorted(s_list,key=lambda x: int(x.name.split("_")[0]))
    count_dict['raw_count'] = pd.DataFrame.from_items([(s.name,s) for s in s_list])
    count_dict['raw_ratio'] = count_dict['raw_count']/count_dict['raw_count'].sum(axis=0)

    # corrected counts -> only calculated for CD45PN or CD3PN or TOT samples
    c_excl = ~meta_table['sort'].isin('45NM,EPCAM,CD137,MAIT'.split(","))
    s_list = []
    for donor in donor_list:
        c_donor = (meta_table[donor_key] == donor) & c_excl
        if np.sum(c_donor)==0:
            continue
        order = str(int(adata.obs['Order'][c_donor][0]))
        age = adata.obs['Age'][c_donor][0]
        stage = adata.obs['Stage'][c_donor][0]

        s_list.append(meta_table[c_donor].groupby(anno_key)['factor'].sum())
        s_list[-1] = s_list[-1].rename("%s_%s_%s_%s"%(order, donor, age, stage))


    s_list = sorted(s_list,key=lambda x: int(x.name.split("_")[0]))
    count_dict['norm_count'] = pd.DataFrame.from_items([(s.name,s) for s in s_list])
    count_dict['norm_ratio'] = count_dict['norm_count']/count_dict['norm_count'].sum(axis=0)
    
    return count_dict


def get_df_for_cell_population(cell_list,cnt_df,save_loc=None,figsize=(6,2)):
    
    from scipy.interpolate import make_interp_spline, BSpline
    from statsmodels.nonparametric.smoothers_lowess import lowess
    
    df = cnt_df[cnt_df.index.isin(cell_list)]
    df = (df/df.sum(axis=0))*100
    df = df.reindex(cell_list)

    plt.figure(figsize=figsize)
    ax = plt.subplot(111)

    for ct in cell_list:
        xlen = len(df.columns)
        x1 = np.arange(xlen)
        y1 = np.array(df[df.index==ct])[0]
        mask = ~np.isnan(y1)
        x1 = x1[mask]
        y1 = y1[mask]
        ax.plot(x1,y1,marker='o', 
                alpha=0.5,linewidth=0.5,label=ct)

    ax.grid(False)
    ax.set_xticks(np.arange(xlen))
    ax.set_xticklabels(df.columns)
    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10) 
            #specify integer or one of preset strings, e.g.
            #tick.label.set_fontsize('x-small') 
            tick.label.set_rotation('vertical')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(bbox_to_anchor=(1.0,0.8))
    if save_loc:
        plt.savefig(save_loc,bbox_inches='tight',format='pdf',dpi=300)
    return df


def get_df_for_cell_population_v2(cell_list,count_dict,style='norm',
                                  save_loc=None,figsize=(6,2),ylim=None,size=20):
    
    from scipy.interpolate import make_interp_spline, BSpline
    from statsmodels.nonparametric.smoothers_lowess import lowess
    
    raw_cnt_df = count_dict['%s_count'%style]
    raw_cnt_df = raw_cnt_df.fillna(0)
    cnt_df = count_dict['%s_ratio'%style]
    
    df = cnt_df[cnt_df.index.isin(cell_list)]
    df = (df/df.sum(axis=0))*100
    df = df.reindex(cell_list)

    plt.figure(figsize=figsize)
    ax = plt.subplot(111)

    counts = np.zeros(len(raw_cnt_df.columns))
    for i, ct in enumerate(cell_list):
        xlen = len(df.columns)
        x1 = np.arange(xlen)
        y1 = np.array(df[df.index==ct])[0]
        cnts = np.array(raw_cnt_df[raw_cnt_df.index==ct])[0]
        mks = np.log10(cnts)*size
        counts+=cnts

        mask = ~np.isnan(y1)
        x1 = x1[mask]
        y1 = y1[mask]
        mks = list(mks[mask])
        
        ax.scatter(x1-0.03*len(cell_list)+0.03*i,y1,s=mks,marker='o',alpha=0.9,label=ct,linewidth=0)  
       # ax.plot(x1,y1,marker='.', markersize=0, 
       #         alpha=0.5,linewidth=0.5,label=ct)


    ax.grid(False)
    ax.set_xticks(np.arange(xlen))
    ax.set_xticklabels([k.split("_")[-2] for k in df.columns])
    #ax.set_xticklabels([' '.join(k.split("_")[-2:]) for k in df.columns])
    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10) 
            #specify integer or one of preset strings, e.g.
            #tick.label.set_fontsize('x-small') 
            tick.label.set_rotation('vertical')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.legend(bbox_to_anchor=(1.3,0.8))
    lgnd = ax.legend(bbox_to_anchor=(1.3,0.8), scatterpoints=1, fontsize=10)

    for i in range(len(lgnd.legendHandles)):
        lgnd.legendHandles[i]._sizes = [10]
    
    if ylim:
        plt.ylim(-5,ylim)
    if save_loc:
        plt.savefig(save_loc,bbox_inches='tight',format='pdf',dpi=300)
    return df,counts

def csr_vappend(a,b):
    """ Takes in 2 csr_matrices and appends the second one to the bottom of the first one. 
    Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
    the first matrix instead of copying it. The data, indices, and indptr still get copied."""

    a.data = np.hstack((a.data,b.data))
    a.indices = np.hstack((a.indices,b.indices))
    a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
    a._shape = (a.shape[0]+b.shape[0],b.shape[1])
    return a

def ravel_index(pos, shape):
    ''' inverse of np.unravel_index'''
    res = 0
    acc = 1
    for pi, si in zip(reversed(pos), reversed(shape)):
        res += pi * acc
        acc *= si
    return res

def flatten(l):
    return [item for sublist in l for item in sublist]

def intersect(*d):
    sets = list(map(set, d))
    result = sets[0]
    for s in sets[1:]:
        result = result.intersection(s)
    return result