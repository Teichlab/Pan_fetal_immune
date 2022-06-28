#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Prepare cell2location reference signatures')
parser.add_argument("h5as", type=str,default=None,help='input h5ad file with reference dataset')
parser.add_argument("output", type=str,default=None,help='folder to write output')
args = parser.parse_args()


import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from scvi.model.base._base_model import BaseModelClass

import scvi
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import issparse
from scvi import _CONSTANTS
from scvi.data._anndata import get_from_registry

import matplotlib as mpl


# parameters
# infile='/nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetTH.exclude_lowQ.addTECs.keep_fetal_TECs.h5ad'
# outdir='SP.minAge14'

infile=args.h5as 
outdir=args.output

genescsv = '../visium.genes.csv' # to filter reference genes by visium, just list of visium geneids

labels_key='anno_c2l'

batch_key='Sample.lanes'
categorical_covariate_keys=['method', 'donor']


#######################
# create output folder
os.mkdir(outdir)
# read data
ref = sc.read(infile)

# for subsetTH that have some NaNs in Sample.lanes
ref.obs['Sample.lanes'] = ref.obs['Sample.lanes'].cat.add_categories("NaN").fillna("NaN")

ref.var['GeneID'] = ref.var['GeneID'].astype('string')
ref.var=ref.var.set_index('GeneID')
print('Raw: cells = '+str(ref.shape[0])+"; genes = " + str(ref.shape[1]))

# Filter genes
genes = pd.read_csv(genescsv)

intersect = np.intersect1d(ref.var.index, genes['ENSEMBL'])
ref = ref[:, intersect].copy()

# filter MT
a = [not(gene.startswith('MT-')) for gene in ref.var['GeneName']]
ref = ref[:,a]

selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
plt.savefig(outdir+'/gene.filter.pdf') 
ref = ref[:, selected].copy()

print('After filtering: cells = '+str(ref.shape[0])+"; genes = " + str(ref.shape[1]))
# train
scvi.data.setup_anndata(adata=ref,
                        batch_key=batch_key,
                        labels_key=labels_key,
                        categorical_covariate_keys=categorical_covariate_keys)
scvi.data.view_anndata_setup(ref)


mod = RegressionModel(ref)

mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True,progress_bar_refresh_rate=0)

# plot ELBO loss history during training, removing first 20 epochs from the plot
fig, ax = plt.subplots()
mod.plot_history(20)
plt.savefig(outdir+'/train.history.pdf')

ref = mod.export_posterior(
    ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
mod.save(outdir+"/rsignatures", overwrite=True)
ref.write(outdir+"/rsignatures/sc.h5ad")

# save signatures
inf_aver = ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = ref.uns['mod']['factor_names']
inf_aver.to_csv(outdir+'/rsignatures/inf_aver.csv')

# ref = sc.read_h5ad(outdir+"/rsignatures/sc.h5ad")
# mod = cell2location.models.RegressionModel.load(outdir+"/rsignatures", ref)

# i didn't manage to plot QC....
def plot_QC1(m,plot,summary_name: str = "means",use_n_obs: int = 1000):
  if use_n_obs is not None:
    ind_x = np.random.choice(m.adata.n_obs, np.min((use_n_obs, m.adata.n_obs)), replace=False)
  else:
    ind_x = None
  m.expected_nb_param = m.module.model.compute_expected(
    m.samples[f"post_sample_{summary_name}"], m.adata, ind_x=ind_x
    )
  x_data = get_from_registry(m.adata, _CONSTANTS.X_KEY)[ind_x, :]
  if issparse(x_data):
    x_data = np.asarray(x_data.toarray())
  
  mu = m.expected_nb_param["mu"]
  data_node = x_data
  plot.hist2d(np.log10(data_node.flatten()+1), np.log10(mu.flatten()+1), bins=50, norm=mpl.colors.LogNorm())
  plot.set_title("Reconstruction accuracy")
  plot.set(xlabel="Data, log10", ylabel="Posterior sample, values, log10")


def plot_QC2(m,plot,summary_name: str = "means",use_n_obs: int = 1000,scale_average_detection: bool = True):
  inf_aver = m.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
  if scale_average_detection and ("detection_y_c" in list(m.samples[f"post_sample_{summary_name}"].keys())):
    inf_aver = inf_aver * m.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
  aver = m._compute_cluster_averages(key="_scvi_labels")
  aver = aver[m.factor_names_]
  plot.hist2d(
    np.log10(aver.values.flatten() + 1),
    np.log10(inf_aver.flatten() + 1),
    bins=50,
    norm=mpl.colors.LogNorm(),)
  plot.set(xlabel="Mean expression for every gene in every cluster", ylabel="Estimated expression for every gene in every cluster")



fig, (ax1,ax2) = plt.subplots(1,2)
plot_QC1(mod,plot=ax1,use_n_obs=10000)
plot_QC2(mod,plot=ax2)
plt.tight_layout()
plt.savefig(outdir+'/train.QC.pdf')

# session information
import sys
for module in sys.modules:
    try:
        print(module,sys.modules[module].__version__)
    except:
        try:
            if  type(modules[module].version) is str:
                print(module,sys.modules[module].version)
            else:
                print(module,sys.modules[module].version())
        except:
            try:
                print(module,sys.modules[module].VERSION)
            except:
                pass
