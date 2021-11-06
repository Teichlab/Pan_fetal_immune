### Prepare anndatas as reference for cell2location analysis ###

import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--subset_organ", type=str,
                    default=None)
parser.add_argument("--min_age", type=int,
                    default=None)
parser.add_argument("--split_stroma", type=bool,
                    default=False)
parser.add_argument("--add_TECs", type=bool,
                    default=False)
parser.add_argument("--keep_fetal_TECs", type=bool,
                    default=False)
args = parser.parse_args()

org = args.subset_organ
age = args.min_age
split_stroma = args.split_stroma
add_tecs = args.add_TECs
keep_fetal_TECs = args.keep_fetal_TECs

def make_c2l_reference(
    ref_adata,
    annotation_obs = 'anno_c2l',
    technical_obs = ["method", "donor"], ## Covariates to regress out
    library_obs = ['Sample.lanes'], ## Covariate for 10x library
    subset_organ=None,
    min_age = None,
    exclude_clusters = None, ## Clusters to exclude
    split_by_organ = None, ## for which clusters should the annotation be split by organ? e.g. just stroma
    min_cells = 20 ## Minimum number of cells to keep a cluster
    ):
    ## Subset by organ
    if subset_organ:
        ref_adata = ref_adata[ref_adata.obs["organ"] == subset_organ].copy()

    ## Subset by age
    if min_age:
        ref_adata = ref_adata[ref_adata.obs["age"] >= min_age].copy()

    ## Exclude low quality clusters
    if exclude_clusters:
        ref_adata = ref_adata[~ref_adata.obs[annotation_obs].isin(exclude_clusters)].copy()

    ## Split selected clusters by organ
    if split_by_organ:
        tosplit_obs = ref_adata.obs[ref_adata.obs[annotation_obs].isin(split_by_organ)].copy()
        organ_anno = tosplit_obs[annotation_obs] + "_" + tosplit_obs['organ']
        ref_adata.obs.loc[tosplit_obs.index, annotation_obs] = organ_anno.values

    ## Remove clusters containing less than n cells
    clus_counts = ref_adata.obs[annotation_obs].value_counts() 
    keep_clus = clus_counts.index[clus_counts >= min_cells] 
    ref_adata = ref_adata[ref_adata.obs[annotation_obs].isin(keep_clus)].copy()

    ## Clean obs
    ref_adata.obs = ref_adata.obs[technical_obs + library_obs + [annotation_obs, "organ", "age"]].copy()

    return(ref_adata)

def save_c2l_reference(params):
    outfile = "PAN.A01.v01.c2l_reference.v2."
    if params["subset_organ"]:
        outfile = outfile + "subset{o}.".format(o=params["subset_organ"])
    if params["split_by_organ"]:
        outfile = outfile + "organ_split_stroma."
    if params["min_age"]:
        outfile = outfile + "minAge{a}.".format(a=params["min_age"])
    if params["exclude_clusters"]:
        outfile = outfile + "exclude_lowQ."
    if params["add_TECs"]:
        outfile = outfile + "addTECs."   
    if params["keep_fetal_TECs"]:
        outfile = outfile + "keep_fetal_TECs."   
    outfile = outdir + outfile + "h5ad"
    
    params.pop('add_TECs')
    params.pop('keep_fetal_TECs')
    ref_adata = make_c2l_reference(adata, **params)
    ref_adata.write_h5ad(outfile)

data_dir="/nfs/team205/ed6/data/Fetal_immune/"
timestamp="20210429"

# Make folder to save outputs
outdir = data_dir + "c2l_scRNA_references/"
if not os.path.exists(outdir):
    os.mkdir(outdir)

### Load reference data
adata = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_raw_count.{t}.h5ad'.format(t=timestamp))
adata.var_names_make_unique()

## Filter maternal contaminants
mat_barcodes = pd.read_csv("~/Pan_fetal_immune/metadata/souporcell_results/maternal_barcodes.csv", index_col=0)
mat_barcodes["x"] = pd.Series([x.split("-1")[0] for x in mat_barcodes['x']])
adata = adata[~adata.obs_names.isin(mat_barcodes["x"])]

## Read annotation groupings
import json
with open('../../metadata/anno_groups.json', 'r') as json_file:
    anno_groups_dict = json.load(json_file)

anno_obs = pd.read_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv", index_col=0)
adata = adata[adata.obs_names.isin(anno_obs.index)].copy()
adata.obs = anno_obs.loc[adata.obs_names].copy()
del anno_obs

## Restrict to organs profiled with visium
adata = adata[adata.obs["organ"].isin(["TH", "SP", "LI"])].copy()
import gc
gc.collect()

## Tissue specific annotations
adata.obs["anno_c2l"] = adata.obs["anno_lvl_2_final_clean"].copy()
# Thymus
TH_annotations = pd.read_csv('/home/jovyan/mount/gdrive/Pan_fetal/annotations/original_files/fetal_thymus_anno.csv')
anno_c2l_th = adata.obs[adata.obs['organ'] == "TH"].copy()
anno_c2l_th['anno_organ'] = np.nan
anno_c2l_th.loc[TH_annotations['index'][TH_annotations['index'].isin(anno_c2l_th.index)],'anno_organ'] = TH_annotations[TH_annotations['index'].isin(anno_c2l_th.index)]['Anno_level_5'].values

rename_cells = anno_c2l_th.index[anno_c2l_th['anno_lvl_2_final_clean'] == "KERATINOCYTE"]
rename_label = anno_c2l_th.anno_organ[anno_c2l_th['anno_lvl_2_final_clean'] == "KERATINOCYTE"].values

adata.obs.loc[rename_cells,"anno_c2l"] = rename_label

print(adata)
## Add additional TECs
if add_tecs:
    th_atlas = '/lustre/scratch117/cellgen/team205/cs42/jovyan_25082021/thymusatlas/HTA07.A01.v02.entire_data_raw_count.h5ad'
    th_atlas_anno = '/lustre/scratch117/cellgen/team205/cs42/jovyan_25082021/thymusatlas/HTA08.v01.A05.Science_human_fig1.h5ad'
    # get TECs
    th_adata = sc.read_h5ad(th_atlas)
    th_adata_anno = sc.read_h5ad(th_atlas_anno, backed='r')
    tec_labels = [x for x in th_adata_anno.obs['Anno_level_fig1'].unique() if "TEC" in x]
    tec_cells = th_adata_anno.obs_names[th_adata_anno.obs['Anno_level_fig1'].isin(tec_labels)]
    if keep_fetal_TECs:
        keep_ages = [x for x in th_adata_anno.obs['Age'].unique() if x.endswith('w')]
        tec_cells = th_adata_anno.obs_names[(th_adata_anno.obs['Anno_level_fig1'].isin(tec_labels)) & (th_adata_anno.obs['Age'].isin(keep_ages))]
    tec_adata = th_adata[tec_cells[~tec_cells.isin(anno_c2l_th.index)]].copy()

    adata.var_names_make_unique()
    tec_adata.var_names_make_unique()

    adata = adata.concatenate(tec_adata, join='outer')
    adata.obs.loc[tec_adata.obs_names + "-1" , 'anno_lvl_2_final_clean'] = th_adata_anno[tec_adata.obs_names].obs["Anno_level_fig1"].values
    adata.obs["anno_c2l"] = adata.obs["anno_lvl_2_final_clean"].copy()
    print(adata)

## Save

lowQ_clusters = ('DOUBLET_IMMUNE_FIBROBLAST',
                     'LOW_Q_INCONSISTENT',
                     'DOUBLET_LYMPHOID_MACROPHAGE',
                     'DOUBLETS_FIBRO_ERY',
                     'DOUBLET_ENDOTHELIUM_ERYTHROCYTE',
                     'DOUBLET_ERY_B',
                     'PLACENTAL_CONTAMINANTS',
                     'DOUBLET')

if split_stroma:
    split_stroma = anno_groups_dict["STROMA"]

params = {
                'subset_organ':org,
                'split_by_organ' : split_stroma,
                'min_age':age,
                'exclude_clusters':lowQ_clusters,
                'add_TECs':add_tecs,
                'keep_fetal_TECs':keep_fetal_TECs
            }

save_c2l_reference(params)
