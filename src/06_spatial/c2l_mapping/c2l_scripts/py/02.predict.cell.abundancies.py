#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Prepare cell2location reference signatures')
parser.add_argument("visium", type=str,default=None,help='visium h5ad file')
parser.add_argument("ref", type=str,default=None,help='reference signatures in csv dormat')
parser.add_argument("output", type=str,default=None,help='folder to write output')
#args = parser.parse_args(['/nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.spleen.h5ad','../ref/subsetSP/rsignatures/inf_aver.csv','subsetSP'])
args = parser.parse_args()


import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from scvi.model.base._base_model import BaseModelClass

import pyro
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

#parameters
batch_key='sample'
alpha=20
max_epochs=50000
# init
vis = sc.read(args.visium)
inf_aver = pd.read_csv(args.ref,index_col=0)
os.mkdir(args.output)

################

intersect = np.intersect1d(vis.var_names, inf_aver.index)

vis = vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=vis,batch_key=batch_key)
scvi.data.view_anndata_setup(vis)

# run on Farm, GPU!
mod = cell2location.models.Cell2location( 
    vis, cell_state_df=inf_aver,
    amortised=False,
    N_cells_per_location=30, 
    detection_alpha=alpha
)

mod.train(max_epochs=max_epochs, 
          batch_size=None,
          train_size=1,
          plan_kwargs={'optim': pyro.optim.Adam(optim_args={'lr': 0.002})},
          use_gpu=True,
          progress_bar_refresh_rate=0)

# plot ELBO loss history during training, removing first 100 epochs from the plot
vis = mod.export_posterior(
    vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

fig, ax = plt.subplots()
mod.plot_history(1000)
fig.legend(labels=['full data training']);
fig.savefig(args.output+'/train.history.pdf') 


# Save model
mod.save(args.output+"/predmodel", overwrite=True)
vis.write(args.output+"/predmodel/sp.h5ad")

mod.plot_QC()
plt.savefig(args.output+'/predict.QC.pdf')


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



