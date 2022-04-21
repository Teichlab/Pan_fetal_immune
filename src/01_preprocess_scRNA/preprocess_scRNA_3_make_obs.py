#!/usr/bin/env python
# coding: utf-8

# ## Fetal Immune Atlas - save labels in anndata
# 
# In this script I load the raw merged anndata file and add the labels unified in `20201230_UniformCellLabels.ipynb`, fixing some inconsistencies on cell names between organs.

import os,sys
import numpy as np 
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata
import scipy
parser = argparse.ArgumentParser()
parser.add_argument("timestamp", help="data timestamp to use")
parser.add_argument("--indir", 
                    default="/nfs/team205/ed6/data/Fetal_immune/",
                    help="folder containing anndata obj")
parser.add_argument("--annot_dir", 
                    default='/home/jovyan/mount/gdrive/Pan_fetal/annotations/',
                    help="folder containing unified annotations")
parser.add_argument("--metadata_path", 
                    default="/home/jovyan/mount/gdrive/Pan_fetal/annotations/manifest_clean_120121.csv",
                    help="folder containing unified metadata")
args = parser.parse_args()

data_dir = args.indir
timestamp = args.timestamp
annot_dir = args.annot_dir
metadata_loc = args.metadata_path

# ### Load merged dataset 
merged_raw = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.h5ad'.format(t=timestamp))


# ### Add cell type labels
# 
# Made uniform in `notebooks/20201230_UniformCellLabels.ipynb`

annot_df = pd.read_csv(annot_dir + "uniform_labels_full.csv", index_col=0)

# Fix names to make them uniform w dataset
def _translate_obs_names(x, organ):
    if organ in ["sp", 'bm']:
        if "FCA" in x:
            obs_name = x.split("_")[3]+ '-' +x.split("_")[5]
        else:
            obs_name = x
    elif organ in ["ki"]:
        if "FCA" in x:
            obs_name = x.split("_")[0] + "-" + x.split("_")[-1].split('-')[0]
        else:
            obs_name = x
    elif organ in ["li", 'ys']:
        obs_name = x.split("_")[3]+'-'+x.split("_")[4]
    elif organ in ["sk"]:
        obs_name = x.split("-")[2]+'-'+x.split("-")[0]
    elif organ in ["gu"]:
        obs_name = "-".join(x.split("-")[:2]) + "_" + x.split("-")[2]
    else:
        obs_name = x
    return(obs_name)

annot_df.index = ["GEX".join(x.split("prime")) for x in annot_df.index]
new_name = [_translate_obs_names(annot_df.index[i],annot_df.organ[i]) for i in range(annot_df.shape[0])]
annot_df["old_name"] = annot_df.index
annot_df.index = new_name

## Subset to cells in the adata
annot_df = annot_df.loc[merged_raw.obs_names[merged_raw.obs_names.isin(annot_df.index)]]

new_anno = pd.concat([merged_raw.obs, annot_df[['uniform_label', 'uniform_label_expanded_merged', 'uniform_label_lvl0']]], 1)
merged_raw.obs = new_anno.loc[merged_raw.obs_names]


### Add metadata
metadata = pd.read_csv(metadata_loc, index_col=0)

## Add library prep info
metadata['method'] = [x.split("prime")[0]+"GEX" if "prime" in x else x for x in metadata["Sequencing"]]

## Rename columns as they are in obs
metadata['donor'] = metadata['SAMPLE.NAME']

## Select useful columns
clean_metadata = metadata[["Organ","Sample.lanes", "Sort_id","age", "method", "donor", "sex", "Processing_method", "AnnatomicalPart"]]
clean_metadata["file"] = clean_metadata['Sample.lanes']

## Rename organs
clean_metadata = clean_metadata.rename({"Organ":"organ"}, axis=1)

organ_ids = {
    "skin":"SK",
    "spleen":"SP",
    "yolkSac":"YS",
    "liver":"LI",
    "thymus":"TH",
    "gut":"GU",
    "boneMarrow":"BM",
    "kidney":"KI"
}

clean_metadata["organ"] = [organ_ids[x] for x in clean_metadata["organ"]]

## Change nomenclature for Mesenteric Lymphnode
clean_metadata['organ'] = clean_metadata.organ.astype("str")
clean_metadata.loc[clean_metadata.index[clean_metadata.AnnatomicalPart=="MLN"],'organ'] = "MLN"

## Save sampleID
clean_metadata["Sample"] = clean_metadata[["donor", "organ", "Sort_id", "file"]].agg('_'.join, axis=1)

## Join metadata with anndata.obs
merged_raw.obs["file"] = ["_".join(x.split("_")[3:]).split("_GRCh")[0] if "cellranger" in x else x for x in merged_raw.obs["file"]]
new_obs = merged_raw.obs.reset_index().merge(clean_metadata, on=['file'], how='left', indicator=True)
new_obs = new_obs.set_index("index")

## Check that the merge has worked properly
if not new_obs.shape[0] == merged_raw.obs.shape[0]:
    print("--- WARNING!! The new obs has more rows than the old obs ---")

if np.any(new_obs._merge=="right_only"):
    print("--- WARNING!! Some values are unique to metadata ---")

if not new_obs.index.is_unique:
    print("--- WARNING!! Duplicate indices ---")
    
new_obs = new_obs.drop(["_merge"],1)

## Save
new_obs.to_csv(data_dir + "PAN.A01.v01.entire_data_normalised_log.{t}.full_obs.csv".format(t=timestamp))

