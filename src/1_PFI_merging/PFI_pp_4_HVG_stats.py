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
args = parser.parse_args()

data_dir = args.indir
timestamp = args.timestamp

adata_lognorm = sc.read_h5ad(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.h5ad'.format(t=timestamp))
adata_lognorm.var_names_make_unique()

## Save dispersion estimates in adata_lognorm for HVG selection
sc.pp.highly_variable_genes(adata_lognorm, min_mean=0.001, max_mean=10, subset=False)
adata_lognorm.var.to_csv(data_dir + 'PAN.A01.v01.entire_data_normalised_log.{t}.var.csv'.format(t=timestamp))