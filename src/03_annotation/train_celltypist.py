import scanpy as sc
import celltypist
import pandas as pd
import numpy as np
import os

######-----------------------------------VERSION1-----------------------------------######
adata = sc.read('/nfs/team205/ed6/data/Fetal_immune/cellxgene_h5ad_files/scRNA_data/PAN.A01.v01.raw_count.20210429.PFI.embedding.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
model = celltypist.train(adata, 'anno_lvl_2_final_clean', use_SGD = True, n_jobs = -1, feature_selection = True, top_genes = 300, details = 'stromal and immune populations from the human fetus', url = 'https://celltypist.cog.sanger.ac.uk/models/Pan_Fetal_Suo/v1/Pan_Fetal_Human.pkl', source = 'internal', version = 'v1')
model.write('../model/v1/Pan_Fetal_Human.pkl')
