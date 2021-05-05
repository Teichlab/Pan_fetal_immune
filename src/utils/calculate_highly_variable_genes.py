## Save gene dispersion estimates ##
# used for HVG selection
import scanpy as sc
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("timestamp", help="data timestamp to use")
parser.add_argument("--indir", 
                    default="/nfs/team205/ed6/data/Fetal_immune/",
                    help="folder containing anndata obj")
parser.add_argument("--split_name", 
                    default="",
                    help="ID for data split (e.g. NKT, Progenitors, Stroma...) (default: no split, full atlas)")
args = parser.parse_args()

data_dir = args.indir
timestamp = args.timestamp
spl = args.split_name
if len(spl)==0:
    h5ad_file = data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.h5ad'.format(t=timestamp)
else:
    h5ad_file = data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.{s}.h5ad'.format(t=timestamp, s=spl)

adata_lognorm = sc.read_h5ad(h5ad_file)
adata_lognorm.var_names_make_unique()

## Save dispersion estimates in adata_lognorm for HVG selection
sc.pp.highly_variable_genes(adata_lognorm, min_mean=0.001, max_mean=10, subset=False)
adata_lognorm.var.to_csv(h5ad_file.split(".h5ad")[0] + ".var.csv")