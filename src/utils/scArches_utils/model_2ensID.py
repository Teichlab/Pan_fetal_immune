## Convert var_names in scvi model to EnsemblIDs for sharing
import pandas as pd
import argparse 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("split_name", help="ID for data split (e.g. NKT, Progenitors, Stroma... PFI for full dataset)")
args = parser.parse_args()

split = args.split_name

scvi_outs_dir = "/home/jovyan/mount/gdrive/Pan_fetal/data4gpu_node/"
model_dir = None
timestamp = "20210429"

if model_dir is None:
        model_dir='scvi_' + split + '_model/'
    ## Check that the feature used for reference training match the var_names in query data
var_names_model = pd.read_csv(scvi_outs_dir + model_dir + "var_names.csv", header=None)[0].values

if split == 'PFI':
    adata_ref_var = pd.read_csv(scvi_outs_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.var.csv'.format(t=timestamp, s=split), index_col=0)
else:
    adata_ref_var = pd.read_csv(scvi_outs_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.var.csv'.format(t=timestamp, s=split), index_col=0)
adata_ref_var.iloc[var_names_model]['GeneID'].to_csv(scvi_outs_dir + model_dir + "var_names.csv", header=None, index=False)