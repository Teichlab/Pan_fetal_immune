import scanpy as sc
import pandas as pd
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma...)")
args = parser.parse_args()

split_name = args.split_name

save_path = '/nfs/team205/ed6/data/Fetal_immune/'
prefix = 'PAN.A01.v01.entire_data_normalised_log.wGut.batchCorrected_20210118'

## Load full datasets (w/o gene selection)
print("Started loading...")
if split_name=="PFI":
    merged_raw = sc.read_h5ad(save_path + 'PAN.A01.v01.entire_data_normalised_log.wGut.h5ad')
    ## Fix obs_names
    def _rename_gut_cells(x):
        if "FCA" not in x:
            x = x.split("_")[8].split('-')[1] + "-1"  + "_" + "_".join(x.split("_")[3:6])
        else: 
            x = x.split("_")[7].split('-')[1] + "-1" + "_" + "_".join(x.split("_")[3:5]) 
        return(x)

    obs_names = merged_raw.obs_names.values
    gut_ixs = np.where(merged_raw.obs.organ=="GU")[0]
    for i in gut_ixs:
        obs_names[i] = _rename_gut_cells(obs_names[i])

    merged_raw.obs_names = obs_names
else:
    merged_raw  = sc.read_h5ad(save_path + '{p}.{s}.h5ad'.format(p=prefix, s=split_name))
print("Finished loading!")
# merged_raw.obs['batch'] = [x+y for x,y in zip(merged_raw.obs['organ'],merged_raw.obs['method'])]
# merged_raw.obs['bbk'] = [x+y for x,y in zip(merged_raw.obs['donor'],merged_raw.obs['method'])]

## Load csv file containing split info
split_info_dir = '/nfs/team205/ed6/data/Fetal_immune/splits_info/'
adata_obs = pd.read_csv(split_info_dir + "{p}.split_info.{s}.csv".format(p=prefix, s=split_name), index_col=0)

merged_raw = merged_raw[adata_obs.index]
## Add info on propagated label 4 ridge regression
merged_raw.obs["uniform_label_propagated"] = adata_obs["uniform_label_propagated"]

## Save split anndata
for split_col in adata_obs.columns:
    if "/" in split_col:
        adata_obs[split_col] = ["_".join(x.split("/")) for x in adata_obs[split_col]]
        adata_obs.columns = ["_".join(x.split("/")) for x in adata_obs.columns]
        split_col = "_".join(split_col.split("/"))
    s = split_col.lstrip("split_")
    adata_name = save_path + "{}.{}.h5ad".format(prefix, s)
    sdata = merged_raw[adata_obs[split_col]==s]
    print("Saving {} anndata ({} cells, {} organs)".format(s, sdata.obs_names.shape[0], sdata.obs["organ"].unique().shape[0]))
    sdata.write_h5ad(adata_name)