# import common modules
from collections import Counter, defaultdict
import os, sys, glob, re
import scipy
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns

# single-cell analysis modules
import scanpy as sc
import scanpy.external as sce
from geosketch import gs
import scrublet as scr
from bbknn import bbknn

# define folders
cwd = '/mnt/18_Pan_fetal/scjp'
data_folder = os.path.join(cwd,'data')
sys.path.append(cwd)

# internal import

from .genes import cc_genes
from .colors import vega_20, vega_20_scanpy, zeileis_26, godsnot_64
from .markers import find_markers, show_marker_plot
from .model import generate_training_X

from . import network
from . import model
from . import species
from . import markers
from . import utils
from . import genes

# matrix files
matrix_files = open(os.path.join(cwd,'./matrix/h5ad_files.txt')).read()

# import genesets
tf_df = pd.read_csv(os.path.join(data_folder,'D01_Human_TF.txt'),sep='\t')
tf_genes = list(tf_df.Symbol) # tf gene list (from Animal TF db)
cd_genes = pd.read_csv(os.path.join(data_folder,'D02_CD_genes.txt'),sep='\t') # CD gene list (HUGO)

# General single-cell anlaysis procedure

def sc_process(adata,pid = 'fspkuc',n_pcs=50): # simplified scanpy preprocessing
    '''n: normalise
       l: log
       f: filter hvg
       r: remove cc_genes
       s: scale
       p: pca
       k: knn_neighbors
       u: umap
       c: leiden clusering
       '''
    if 'n' in pid:
        sc.pp.normalize_per_cell(adata,counts_per_cell_after=10e4)
    if 'l' in pid:
        sc.pp.log1p(adata)
        adata.raw = adata
        print('adding raw...')
    if 'f' in pid:
        if adata.raw == None:
            adata.raw = adata
            print('adding raw...')
        sc.pp.filter_genes_dispersion(adata)
    if 'r' in pid:
        adata = remove_geneset(adata,cc_genes)
        print('removing cc_genes...')
    if 's' in pid:
        sc.pp.scale(adata,max_value=10)
    if 'p' in pid:
        sc.pp.pca(adata)
    if 'k' in pid:
        sc.pp.neighbors(adata,n_pcs=n_pcs)
    if 'u' in pid:
        sc.tl.umap(adata)
    if 'c' in pid:
        sc.tl.leiden(adata)
    return adata
    
def read_process(adata,version,
                 species = 'human',
                 sample = None,
                 define_var = True,
                 call_doublet = True,
                 write = True,
                 min_n_counts = 1000,
                 min_n_genes = 500,
                 max_n_genes = 7000,
                 max_p_mito = 0.5
                ):
    if sample:
        adata.obs['Sample'] = sample
    if define_var:
        adata.var['GeneName'] = list(adata.var.gene_ids.index)
        adata.var['EnsemblID'] = list(adata.var.gene_ids)
    adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
    adata.obs['n_genes'] = np.sum(adata.X>0,axis=1)
    
    print('calculating mito... as species = {}'.format(species))
    if species=='mouse':
        mito_genes = adata.var_names.str.startswith('mt-')
        adata.obs['mito'] = (np.sum(adata.X[:, mito_genes],axis=1).A1) / ((np.sum(adata.X,axis=1).A1)+1)
    elif species=='human':
        mito_genes = adata.var_names.str.startswith('MT-')
        adata.obs['mito'] = (np.sum(adata.X[:, mito_genes],axis=1).A1) / ((np.sum(adata.X,axis=1).A1)+1)
    else:
        print("check_species: weird mito gene names")
        raise SystemError
        
    print('filtering cells... higher than {} counts, more than {} and less than {} genes, less than {} p_mito...'.format(min_n_counts, min_n_genes, max_n_genes, max_p_mito))
    # filter cells
    clist = []
    clist.append(np.array(adata.obs['n_counts'] > min_n_counts))
    clist.append(np.array(adata.obs['n_genes'] > min_n_genes))
    clist.append(np.array(adata.obs['n_genes'] < max_n_genes))
    clist.append(np.array(adata.obs['mito'] < max_p_mito))
    
    c = np.column_stack(clist).all(axis=1)
    adata = adata[c].copy()
    
    if call_doublet:
        print('calling doublets using scrublet...')
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets
    
    if write:
        print('writing output into write/%s%s_filtered.h5ad ...'%(version,sample))
        sc.write('%s%s_filtered'%(version,sample),adata)
    return adata

def write_notebook(name1,name2):
    os.system(f'jupyter nbconvert {name1} --to notebook --ClearOutputPreprocessor.enabled=True --output {name2}')
        
def sort_var_names_based_on_GeneID(adata):
    return adata[:,np.argsort(adata.var.GeneID)].copy()

def combine_batch(adata,key1,key2,new_key = 'batch'):
    print('storing new batch into '+new_key)
    adata.obs[new_key] = ['{}_{}'.format(k1,k2) for k1,k2 in zip(adata.obs[key1],adata.obs[key2])]

def doublet(adata, key='Sample'):
    '''detecting doublet using scrublet per key'''
    doublet = []
    for filename in set(adata.obs[key]):
        print(filename)
        sdata = adata[adata.obs[key] == filename].copy()
        scrub = scr.Scrublet(sdata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        doublet.extend([(x,y,z) for x,y,z in zip(sdata.obs_names,doublet_scores,predicted_doublets)])
    doublet_score = {x:y for (x,y,z) in doublet}
    doublet_predict = {x:z for (x,y,z) in doublet}
    adata.obs['doublet_score'] = [doublet_score[obs_name] for obs_name in list(adata.obs_names)]
    adata.obs['doublet_predict'] = [doublet_predict[obs_name] for obs_name in list(adata.obs_names)]

def get_sketch(adata,key,folds=10,how='pd',min_num_per_key=500,start='filter',raw=True):
    '''geometric sketching based on diffusion map and pca
    folds: folds to subsample
    min_num_per_key: minimun number to sample'''
    sketch_index = []
    for smp in set(adata.obs[key]):
        print(smp)
        c = adata.obs[key] == smp
        print('from:',sum(c))
        
        if start=='filter':
            sdata = get_subset(adata,c,raw=raw)
        else:        
            sdata = adata[c]
            sc.pp.filter_genes_dispersion(sdata)
            sc.pp.pca(sdata)
        
        if 'd' in how:
            sc.pp.neighbors(sdata)
            sc.tl.diffmap(sdata)

        N = np.max([np.int(np.sum(c)/folds),np.min([min_num_per_key,np.sum(c)])])
        print('to select:',N)
        if how =='pd':
            set1 = set(sdata.obs_names[gs(sdata.obsm['X_diffmap'],N,replace=False)])
            set2 = set(sdata.obs_names[gs(sdata.obsm['X_pca'][:,:50],N,replace=False)])
            sketch_index.extend(list(set1.union(set2)))
        elif how =='p':
            set2 = set(sdata.obs_names[gs(sdata.obsm['X_pca'][:,:50],N,replace=False)])
            sketch_index.extend(list(set2))
        elif how =='d':
            set1 = set(sdata.obs_names[gs(sdata.obsm['X_diffmap'][:,:20],N,replace=False)])
            sketch_index.extend(list(set1))
        else:
            raise SystemError
        print('length of sketch:',len(sketch_index))
    return(sketch_index)
        
def bbknn_umap(adata,batch_key,n_pcs,cluster=False,n_neighbors=3):
    bbknn(adata,batch_key=batch_key,n_pcs=n_pcs,approx=False,neighbors_within_batch=n_neighbors)
    if cluster:
        sc.tl.leiden(adata)
    sc.tl.umap(adata)
    
def umap(adata,name=None):
    sc.tl.umap(adata)
    if name:
        adata.obsm['X_umap_'+name] = adata.obsm['X_umap'].copy()
        
def umap_show(adata,feature,feature_name= None):
    if feature_name:
        adata.obs[feature_name] = feature
        sc.pl.umap(adata,color=feature_name,color_map='OrRd')
    else:
        adata.obs['show'] = feature
        sc.pl.umap(adata,color='show',color_map='OrRd')

# Clustering
        
def leiden_res(adata,res,show=False):
    print('calculating leiden at res {0:.2f}...'.format(res))
    sc.tl.leiden(adata,resolution=res)
    print('copying into obs.leiden_{0:.2f}...'.format(res))
    adata.obs['leiden_{0:.2f}'.format(res)] = adata.obs['leiden'].copy()
    if show:
        sc.pl.umap(adata,color='leiden_{0:.2f}'.format(res))
        
def subcluster(adata,obs_label,cl_label,new_label,res=0.1):
    '''
    take specific cluster from adata and split that into smaller cluster
    adata: AnnData object
    obs_label: obs label. eg. 'leiden' or 'celltype'
    cl_label: cluster name. eg. '1' or 'macrophage'
    new_label: name to store updated label
    '''
    subset = adata[adata.obs[obs_label]==cl_label].copy()
    sc.tl.leiden(subset,resolution=res)
    update_dict = {obs_name:cl_label+'_'+new_cl for obs_name,new_cl in zip(subset.obs_names,subset.obs['leiden'])}
    adata.obs[new_label] = [old if obs_name not in update_dict else update_dict[obs_name] for obs_name, old in zip(adata.obs_names,adata.obs[obs_label])]
        
def remove_geneset(adata,geneset):
    adata = adata[:,~adata.var_names.isin(list(geneset))].copy()
    return adata

def is_cycling(adata,cc_genes=cc_genes,cut_off=0.4):
    X = np.mean(adata.raw.X[:,adata.raw.var_names.isin(cc_genes)],axis=1)
    plt.hist(X)
    adata.obs['Cycle_score'] = X
    adata.obs['isCycle'] = X>cut_off
    
def get_subset(idata, select, cc_genes=cc_genes, log=False,raw=True):
    if raw:
        adata = sc.AnnData(idata[select].raw.X)
        adata.var = idata.raw.var
    else:
        adata = sc.AnnData(idata[select].X)
        adata.var = idata.var
    adata.obs = idata.obs[select]
    adata.raw = adata.copy()
    #adata.X = scipy.sparse.csr_matrix(np.exp(adata.X.todense())-1)
    sc.pp.filter_genes_dispersion(adata,log=log)
    if log:
        sc.pp.log1p(adata)
    sc.pp.scale(adata,max_value=10)
    if len(cc_genes)>0:
        print('removing cc_genes...')
        adata = remove_geneset(adata,cc_genes)
    sc.pp.pca(adata,n_comps = np.min([50,adata.X.shape[0],adata.X.shape[1]]))
    return adata

def get_raw(idata):
    adata = sc.AnnData(idata.raw.X)
    adata.var = idata.raw.var
    adata.obs = idata.obs
    adata.raw = adata.copy()

    return adata

def get_raw_process(idata, cc_genes=cc_genes, log=False):
    adata = sc.AnnData(idata.raw.X)
    adata.var = idata.raw.var
    adata.obs = idata.obs
    adata.raw = adata.copy()
    #adata.X = scipy.sparse.csr_matrix(np.exp(adata.X.todense())-1)
    sc.pp.filter_genes_dispersion(adata,log=log)
    if log:
        sc.pp.log1p(adata)
    sc.pp.scale(adata,max_value=10)
    if len(cc_genes)>0:
        print('removing cc_genes...')
        adata = remove_geneset(adata,cc_genes)
    sc.pp.pca(adata,n_comps = np.min([50,adata.X.shape[0],adata.X.shape[1]]))
    return adata

def output_matrix_Seurat(adata,version,name,use_raw=False):
    from scipy.io import mmwrite
    
    if use_raw:
        X = adata.raw.X
        mmwrite(version+name+'.mtx',X)
        adata.obs.to_csv(version+name+'.meta.csv')
        adata.raw.var.to_csv(version+name+'.var.csv')
    else:
        X = adata.X
        mmwrite(version+name+'.mtx',X)
        adata.obs.to_csv(version+name+'.meta.csv')
        adata.var.to_csv(version+name+'.var.csv')

def us(adata,gene,groups=None, show=False, exclude =None,figsize=None,**kwargs):
    from matplotlib import rcParams
    if figsize:
        rcParams['figure.figsize'] = figsize
    if ',' in gene:
        gene = gene.split(',')
    if groups:
        sc.pl.umap(adata,color=gene,color_map='OrRd',groups=groups, show=show, **kwargs)
    else:
        if exclude: # list to not show
            groups = [x for x in set(adata.obs[gene]) if x not in exclude]
            sc.pl.umap(adata,color=gene,color_map='OrRd',groups=groups, show=show, **kwargs)
        else:
            sc.pl.umap(adata,color=gene,color_map='OrRd',show=show, **kwargs)
    rcParams['figure.figsize'] = [5,5]
    
def merge_matrix(ad,obskeys = None,use_raw = False,keep_only_mutual=False):
    '''merge matrix stored in ad
    ad: dictionary of anndata to merge
    obskeys: list to merge within anndata
    use_raw: if True, merge from .raw.X'''
    
    smp_list = list(ad.keys())
    obs_dict = defaultdict(list)
    obs_names = []
    
    for smp in smp_list:
        ad[smp].obs['name'] = smp
    
    if not obskeys:
        obskey_list = []
        obskeys = []
        for sample in smp_list:
            obskey_list.extend(list(ad[sample].obs.columns))
        for (obskey, number) in Counter(obskey_list).items():
            if number == len(smp_list):
                obskeys.append(obskey)
            else:
                if keep_only_mutual:
                    pass
                else:
                    for sample in smp_list:
                        if obskey not in ad[sample].obs.columns:
                            ad[sample].obs[obskey]='n/a'
                    obskeys.append(obskey)
                               
    for sample in smp_list:
        obs_names.extend(list(ad[sample].obs_names))
        for key in obskeys:   
            obs_dict[key].extend(list(ad[sample].obs[key]))
    
    from scipy.sparse import vstack
    if use_raw == True:
        stack = vstack([ad[x].raw.X for x in smp_list]) # stack data
        adata = sc.AnnData(stack, var = ad[smp_list[0]].raw.var)
    else:
        stack = vstack([ad[x].X for x in smp_list]) # stack data
        adata = sc.AnnData(stack, var = ad[smp_list[0]].var)
      
    
    adata.obs_names = obs_names
    print(len(adata))
    for obs_col in obs_dict:
        print(obs_col)
        adata.obs[obs_col] = obs_dict[obs_col]
    return adata

def timestamp():
    from datetime import datetime
    return datetime.now().strftime("%y%m%d%H%M")

def save_html(name,log=False): # export notebook into html file
    time = timestamp()
    print(os.system('jupyter nbconvert --to html %s'%(name)))
    name_key = re.sub('.ipynb$','',name)
    if log:
        print(os.system('mv %s.html %s_%s_%s.html'%(name_key,name_key,time,log)))
    else:
        print(os.system('mv %s.html %s_%s.html'%(name_key,name_key,time)))

def write(adata,version,name):
    '''write adata into [name]'''
    name = version + name
    sc.write(name,adata)
    print("_".join(name.split(".")) + " = '%s'"%name)
    
def save_fig(version,figcount,fig_format='pdf',fig_folder='11_Figs'):
    
    plt.savefig('%s/%s%s.%s'%(fig_folder,version,figcount,fig_format),bbox_inches='tight',format=fig_format,dpi=300)
    print('%s/%s%s.pdf'%(fig_folder,version,figcount))

# batch regression methods

def regress_batch_v2(adata,batch_key,confounder_key):
    '''batch regression tool
    batch_key=list of observation categories to be regressed out
    confounder_key=list of observation categories to be kept
    returns ndata with corrected X'''

    from sklearn.linear_model import Ridge
    
    print('fitting linear model...')
    dummy = pd.get_dummies(adata.obs[batch_key+confounder_key],drop_first=False)
    X_exp = adata.X # scaled data
    if scipy.sparse.issparse(X_exp):
        X_exp = X_exp.todense()
    LR = Ridge(fit_intercept=False,alpha=1.0)
    LR.fit(dummy,X_exp)

    if len(batch_key)>1:
        batch_index = np.logical_or.reduce(np.vstack([dummy.columns.str.startswith(x) for x in batch_key]))
    else:
        batch_index = np.vstack([dummy.columns.str.startswith(x) for x in batch_key])[0]
    
    print('corrcting batch...')
    dm = np.array(dummy)[:,batch_index]
    X_explained = dm.dot(LR.coef_[:,batch_index].T)
    X_remain = X_exp - X_explained
    ndata = sc.AnnData(X_remain)
    ndata.obs = adata.obs
    ndata.var = adata.var
    return ndata, X_explained

def regress_iter(adata,batch_key,confounder_key,bbknn_key,scale=True, approx = True,n_pcs=50):
    if scale == True:
        print('scaling data...')
        sc.pp.scale(adata,max_value=10)
    ndata, X_explained = regress_batch_v2(adata,batch_key=batch_key,confounder_key=confounder_key)
    print('running pca...')
    sc.pp.pca(ndata)
    print('running bbknn...')
    bbknn(ndata, batch_key = bbknn_key,n_pcs=n_pcs, approx=approx)
    return ndata #, X_explained

def run_pca_bbknn_umap(ad,level_key,bbknn_key,marker_dict,
                       resolution=0.02,start = 'leiden',select=False,
                      thres=0.95,min_drop_cut=0.5,show=True, how='almost',
                      min_cluster_num = 200):
    '''
    run pca, bbknn, umap and clustering to find good low-rescluster with markers
    how = [almost, any, all] for marker_found function
    '''
    adata = ad[level_key]
    if start == 'pca':
        sc.pp.pca(adata)
        bbknn(adata,batch_key=bbknn_key,approx=False)
        sc.tl.umap(adata)
    if start in 'leiden,pca'.split(','):
        sc.tl.leiden(adata,resolution=resolution)
        adata.obs[level_key] = list(adata.obs['leiden'])
    if start in 'leiden,pca,mks'.split(','):
        if show:
            sc.pl.umap(adata,color='leiden')
        if len(set(adata.obs['leiden']))<2:
            print('clustering not enough')
            return False,None
        elif np.min([x for x in Counter(adata.obs['leiden']).values()]) < min_cluster_num:
            print('clustering resolution too high')
            return False,None
        else:
            if select:
                ndata = generate_training_X(adata,'leiden',select_num=select)
            else:
                ndata = adata
            mks = find_markers(ndata,'leiden',thres=thres,min_drop_cut=min_drop_cut)
            if show:
                show_marker_plot(ndata,'leiden',mks,toshow=5)
            go = marker_found(mks,how=how)
            if go:
                commit_level(adata,level_key,mks,marker_dict)
            else:
                print('marker not found')
            return go, mks
    else:
        'start accepts either pca or leiden'
        raise SystemError

def marker_found(mks,how='any'):
    if how=='any':
        c1 = len([keys for keys,values in mks.items() if len(values)>0]) > 0 
    elif how == 'some':
        c0 = len([keys for keys,values in mks.items() if len(values)>0]) >= 3
        c2 = len([keys for keys,values in mks.items() if len(values)>0]) >= (len(mks.keys())-1)
        c1 = c0|c2
    elif how=='all':
        c1 = len([keys for keys,values in mks.items() if len(values)>0]) == len(mks.keys())
    elif how=='almost':
        c1 = len([keys for keys,values in mks.items() if len(values)>0]) >= (len(mks.keys())-1)
    else:
        print('Error: print how not in any, all, alomst')
        raise SystemExit
    return c1

def commit_level(adata,level_key,mks,marker_dict):
    for leiden_clst in set(adata.obs[level_key]):
        to_merge = np.array(adata.obs[level_key].copy(),dtype=object)
        if len(mks[leiden_clst]) >0:
            final_key = level_key+"_"+leiden_clst
            marker_dict[final_key] = mks[leiden_clst]
        else:
            to_merge[adata.obs[level_key]==leiden_clst]='M'
        adata.obs[level_key] = to_merge
                
def expand_level_copy(ad,level_key):
    adata = ad[level_key]
    for leiden_clst in set(adata.obs[level_key]):
        final_key = level_key+"_"+leiden_clst
        print(final_key)
        ad[final_key] = adata[adata.obs[level_key]==leiden_clst].copy()
        
def summary(ad):
    anno = np.zeros(len(ad['0']),dtype=object)
    final_clusters = []
    for k in ad.keys():
        if np.sum([x.startswith(k) for x in ad.keys()]) ==1:
            final_clusters.append(k)
    for k in final_clusters:
        anno[ad['0'].obs_names.isin(ad[k].obs_names)] =k
    return anno
        
def walk_cluster(ad,marker_dict,tried,bbknn_key,
                 leiden_walk=[0.02,0.05], thres=0.95, min_drop_cut=0.5,
                 select=False, show=False, how='almost', 
                 final_limit_num=8, min_num_split=500):
    go = False
    processed = set(['_'.join(x.split('_')[:-1]) for x in marker_dict.keys()])
    to_process = [level_key for level_key in list(ad.keys()) if level_key not in processed.union(set(tried))]

    print(to_process)
    for level_key in to_process:
        print(level_key)
        if len(level_key.split("_")) > final_limit_num:
            print('level too deep')
            continue
        if len(ad[level_key])<min_num_split:
            print('subset too small')
            continue
        for resolution in leiden_walk:
            print(resolution)
            result = run_pca_bbknn_umap(ad,level_key,bbknn_key,marker_dict,
                                             start='leiden',resolution=resolution,
                                            thres=thres,min_drop_cut=min_drop_cut,
                                            select=select,show=show,how=how)
            if result[0]:
                print('marker found at '+str(resolution))
                go = True
                expand_level_copy(ad,level_key)
                break
        tried.append(level_key)
    return(go)

def run_harmony(sdata,vars_use,n_pcs=30):
    meta_data = sdata.obs
    data_mat = sdata.obsm['X_pca']

    import harmonypy as hm
    ho = hm.run_harmony(data_mat,meta_data,vars_use)
    
    sdata.obsm['X_pca_before'] = sdata.obsm['X_pca']
    sdata.obsm['X_pca'] = ho.Z_corr.T
    
    print('calculating neighbors...')
    sc.pp.neighbors(sdata,n_pcs=n_pcs)
    print('calculating umap...')
    sc.tl.umap(sdata)

# easy annotation class
        
class annotater():
    '''
    create de novo annotation onto adata
    '''
    
    def __init__(self,adata,new_label_name,old_label=None):
    
        if old_label:
            adata.obs[new_label_name] = adata.obs[old_label]
        else:
            adata.obs[new_label_name] = 'unknown'
        arr = np.array(adata.obs[new_label_name],dtype=object)
        self.new_label = arr
        self.new_label_name = new_label_name
        
    def update(self,adata,obskey,select,label_name,unknown=False):
        if type(select)==str:
            if ',' in select:
                label_condition = adata.obs[obskey].isin(select.split(','))
            else:
                label_condition = adata.obs[obskey]==select
        else:
            label_condition = adata.obs[obskey]==select
        if unknown:
            print('updating only unknown values...')
            label_condition = label_condition & (adata.obs[self.new_label_name]=='unknown')
        self.new_label[label_condition] = label_name
        adata.obs[self.new_label_name] = self.new_label
    
    def update_condi(self,adata,label_condition,label_name):

        self.new_label[label_condition] = label_name
        adata.obs[self.new_label_name] = self.new_label

def read_files(filename, sample):
    '''import 10X data, based on filename (path to file) and sample ID (assigned as unique ID)'''

    path = '%s/'%(filename)
    adata = sc.read(path+'matrix.mtx',cache=True).transpose()
  
    try: 
        adata.var_names = np.genfromtxt(path + 'genes.tsv',dtype=str)[:,1]
        adata.var['GeneName'] = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 1]
        adata.var['GeneID'] = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 0]
    except: 
        adata.var_names = np.genfromtxt(path + 'features.tsv.gz',dtype=str)[:,1]
        adata.var['GeneName'] = np.genfromtxt(path + 'features.tsv.gz', dtype=str)[:, 1]
        adata.var['GeneID'] = np.genfromtxt(path + 'features.tsv.gz', dtype=str)[:, 0]
    adata.obs_names = np.genfromtxt(path + 'barcodes.tsv',dtype=str)
    adata.obs_names = [filename+"-"+x.strip("-1") for x in adata.obs_names]
    adata.obs['Sample'] = sample
    
    # caculate n_counts / n_genes per cell
    adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
    adata.obs['n_genes'] = np.sum(adata.X>0,axis=1)
    
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['mito'] = (np.sum(adata.X[:, mito_genes],axis=1).A1) / (np.sum(adata.X,axis=1).A1)
    
    # filter cells
    clist = []
    clist.append(np.array(adata.obs['n_counts'] > 1000))
    clist.append(np.array(adata.obs['n_genes'] > 500))
    clist.append(np.array(adata.obs['n_genes'] < 7000))
    clist.append(np.array(adata.obs['mito'] < 0.5))

    c = np.column_stack(clist).all(axis=1)
    adata = adata[c].copy()

    sc.write('%s%s_filtered'%(version,sample),adata)
    
    return adata

def read_files_multi(file_lists, n_pool = 10):
    '''file_lists: list of tuple(filename, sampleID)'''
    
    import multiprocessing as mp

    with mp.Manager() as manager:
        with manager.Pool(n_pool) as pool:
            pool.starmap(read_files, file_lists)

# Remove doublets

def final_doublets(adata,doublet_key = 'doublet_final',fracDoublet = 0.1,leiden_key = 'leiden'):
    doublet_annotator = annotater(adata, doublet_key)

    for cl in set(adata.obs[leiden_key]):
        isCluster = adata.obs[leiden_key]==cl
        nTotal = np.sum(isCluster)
        try: nDoublet = np.sum(adata.obs['predicted_doublets'][isCluster])
        except: nDoublet = np.sum(adata.obs['predicted_doublets'][isCluster]=='True')
        print(nDoublet,nTotal)
        if nDoublet > fracDoublet*nTotal:
            doublet_annotator.update(adata,leiden_key,cl,'doublet')

# Cell composition analysis

def get_crosstab(adata,tab1,tab2):
    df = pd.crosstab(adata.obs[tab1],adata.obs[tab2])
    df_norm1 = df.div(df.sum(axis=0), axis=1)*100
    df_norm2 = df_norm1.div(df_norm1.sum(axis=1), axis=0)*100
    return df_norm1,df_norm2
            
# linear regression for diff_exp
        
class linear_regression():
    
    def __init__(self,adata):
        self.adata = adata
        self.df = adata.obs
        self.X_exp = adata.X
        self.LR_dict = {}
        
    def ridge(self,keys):
        from sklearn.linear_model import Ridge
        cat = self.df[keys]
        LR = Ridge()
        dummy = pd.get_dummies(cat,drop_first=False)
        LR.fit(dummy, self.X_exp)
        print('fitting linear model...')
        params = list(dummy.columns)
        ct_dict = {'LR': LR, 'params': params}
        self.LR_dict['ridge'] = ct_dict
        
    def lasso(self,keys):
        from sklearn.linear_model import Lasso
        cat = self.df[keys]
        LR = Lasso()
        dummy = pd.get_dummies(cat,drop_first=False)
        LR.fit(dummy, self.X_exp)
        print('fitting linear model...')
        params = list(dummy.columns)
        ct_dict = {'LR': LR, 'params': params}
        self.LR_dict['lasso'] = ct_dict
        
    def celltype_key(self, anno_key, celltype, others=['organ','method']):
        from sklearn.linear_model import Ridge
        
        self.df['is_%s'%celltype] = [str(x) for x in self.df[anno_key]==celltype] # cell type term
        cat = self.df[['is_%s'%celltype] + others] # add organ term and others
        LR = Ridge()
        print('fitting linear model...')
        dummy = pd.get_dummies(cat, drop_first=False)
        LR.fit(dummy, self.X_exp)
        params = list(dummy.columns)
        ct_dict = {'LR': LR, 'params': params}
        self.LR_dict[celltype] = ct_dict
  
        
    def celltype_organ(self, anno_key, celltype, others=['organ','method']):
        from sklearn.linear_model import Ridge
        
        self.df['is_%s'%celltype] = [str(x) for x in self.df[anno_key]==celltype] # cell type term
        self.df['%s_organ'%celltype] = [x+'_'+y for x,y in zip(self.df['is_%s'%celltype],self.df['organ'])] # interaction term
        cat = self.df[['is_%s'%celltype, '%s_organ'%celltype] + others] # add organ term and others
        LR = Ridge()
        dummy = pd.get_dummies(cat, drop_first=False)
        LR.fit(dummy, self.X_exp)
        params = list(dummy.columns)
        ct_dict = {'LR': LR, 'params': params}
        self.LR_dict[celltype] = ct_dict
        
    def celltype_organ_2(self, anno_key, celltype, others=['method']):
        from sklearn.linear_model import Ridge
        
        self.df['is_%s'%celltype] = [str(x) for x in self.df[anno_key]==celltype] # cell type term
        self.df['%s_organ'%celltype] = [x+'_'+y for x,y in zip(self.df['is_%s'%celltype],self.df['organ'])] # interaction term
        cat = self.df[['%s_organ'%celltype] + others] # add organ term and others
        LR = Ridge()
        dummy = pd.get_dummies(cat, drop_first=False)
        LR.fit(dummy, self.X_exp)
        params = list(dummy.columns)
        ct_dict = {'LR': LR, 'params': params}
        self.LR_dict[celltype] = ct_dict
    
    def show_param_genes(self, celltype, param, toshow=20, output=False):
        coef = self.LR_dict[celltype]['LR'].coef_[:,self.LR_dict[celltype]['params'].index(param)]
        order = np.argsort(-coef)
        names = self.adata.var_names[order][:toshow]
        values = coef[order][:toshow]
        print(['%s:%.2f'%(x,y) for x,y in zip(names,values)])
        if output:
            return names, values
        
    def get_param_values(self, celltype, param):
        coef = self.LR_dict[celltype]['LR'].coef_[:,self.LR_dict[celltype]['params'].index(param)]
        return coef
    
    def param_summary(self, celltype, gene, show=True,**kwargs):
        LR = self.LR_dict[celltype]['LR']
        params = self.LR_dict[celltype]['params']

        gidx = self.adata.var_names==gene

        cf = LR.coef_[gidx][0]
        cf_idx = np.argsort(-cf)
        print('\n'.join(['%.2f : %s'%(a,b) for a,b in zip(cf[cf_idx],np.array(params)[cf_idx])]))

    def violin_plot(self, celltype, gene, key='organ',show=True,**kwargs):
        LR = self.LR_dict[celltype]['LR']
        params = self.LR_dict[celltype]['params']
        
        Exp = self.adata.raw.X
        exp = Exp[:,self.adata.raw.var_names==gene].todense().A1

        self.df[gene+'_exp'] = exp

        fig = plt.figure(figsize=(8,2))
        ax = plt.subplot(111)
        sns.violinplot(x='is_%s'%(celltype),y=gene+'_exp',
                       hue=key,data = self.df,**kwargs,
                       scale='width',linewidth=0,inner=None,rasterized=True,cut=0,ax=ax)
        plt.xticks(rotation=0)
        plt.grid(False)
        ax.legend(bbox_to_anchor=(1.2, 1.05))
        if show:
            plt.show()

        gidx = self.adata.var_names==gene

        cf = LR.coef_[gidx][0]
        cf_idx = np.argsort(-cf)
        print('\n'.join(['%.2f : %s'%(a,b) for a,b in zip(cf[cf_idx],np.array(params)[cf_idx])]))
