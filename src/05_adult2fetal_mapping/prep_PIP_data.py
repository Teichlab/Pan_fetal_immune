### Prep pan-immune datasets for mapping ###
import sys,os
import scvi
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

pi_adata = sc.read_h5ad('/nfs/team205/cx1/Celltypist/revision_science/data/PIP_global_object_raw_count.h5ad')
pi_adata.var_names_make_unique()

## Convert gene names to EnsemblIDs as `var_names` (matching with fetal obj)
adata_fetal = sc.read_h5ad('/nfs/team205/ed6/data/Fetal_immune/PAN.A01.v01.entire_data_normalised_log.20210429.h5ad', backed='r')
adata_fetal.var_names_make_unique()
common_genes = pi_adata.var_names[pi_adata.var_names.isin(adata_fetal.var_names)]
print(f'found {len(common_genes)}/{pi_adata.n_vars} common genes')
pi_adata = pi_adata[:,common_genes].copy()
pi_adata.var = adata_fetal.var.loc[common_genes]
pi_adata.var_names = pi_adata.var['GeneID'].values.copy()

## Match column name for 10x method and donor
tr_obs = {'donor_id':'donor', 'chem':"method"}
pi_adata.obs.columns = [tr_obs[x] if x in tr_obs.keys() else x for x in pi_adata.obs.columns]
tr_method = {'5v1':'5GEX', '5v2':'5GEX', '3':'3GEX'}
pi_adata.obs['method'] = [tr_method[x] for x in pi_adata.obs['method']]

### Subset to lymphoid/myeloid cells
lym_anno = ['Tnaive/CM_CD4', 'Memory B cells', 'Trm_gut_CD8', 'NK_CD16+', 'Teffector/EM_CD4',
                   'Trm_Th1/Th17', 'Tfh',
                  'Tem/emra_CD8', 'Naive B cells', 'Trm/em_CD8', 'Tregs',
                  'NK_CD56bright_CD16-', 'Tnaive/CM_CD8', 'Trm_Tgd',
                  'Plasma cells', 'MAIT', 'Tgd_CRTAM+',
                  'Tnaive/CM_CD4_activated', 
                  'Cycling T&NK', 'ILC3', 'ABCs', 'Cycling',
                  'GC_B (I)', 'GC_B (II)',
                  'Pre-B', 'Pro-B', 'Plasmablasts']

mye_anno = ['Classical monocytes','Alveolar macrophages','Mast cells',
            'Nonclassical monocytes', 'Intermediate macrophages','Erythrophagocytic macrophages', 'Progenitor','DC2', 'pDC', 'Intestinal macrophages','DC1', 'Megakaryocytes', 'migDC']

lymph_pi_adata = pi_adata[pi_adata.obs["Category"].isin(lym_anno)].copy()
mye_pi_adata = pi_adata[pi_adata.obs["Category"].isin(mye_anno)].copy()


lymph_pi_adata.write_h5ad("/nfs/team205/ed6/data/Fetal_immune/panimmune_full_LYMPHOID_query.h5ad")
lymph_pi_adata.write_h5ad("/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/panimmune_full_LYMPHOID_query.h5ad")

mye_pi_adata.write_h5ad("/nfs/team205/ed6/data/Fetal_immune/panimmune_full_MYELOID_query.h5ad")
mye_pi_adata.write_h5ad("/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/panimmune_full_MYELOID_query.h5ad")