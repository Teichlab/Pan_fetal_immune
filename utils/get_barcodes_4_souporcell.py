## Save barcodes for souporcell
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np

indir = '/nfs/team205/ed6/data/Fetal_immune/cellbender_raw/'
outdir = '/nfs/team205/ed6/data/Fetal_immune/barcodes4souporcell/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

## Read processed cellranger outputs
# (outputs of PFI_pp_1_read_cellbender.py)
for filename in os.listdir(indir):
    sample = filename.split('_filtered.h5ad')[0]
    sample_outdir = outdir + sample + "/"
    if not os.path.exists(sample_outdir):
        os.mkdir(sample_outdir)
    adata = sc.read_h5ad(indir + filename)
    barcodes = [x.split(sample + "-")[1] + "-1" for x in adata.obs_names]
    f = open( sample_outdir + 'updated_barcodes.tsv','w')
    f.writelines([line+'\n' for line in barcodes])
    f.close()
    
### Save donor 2 sample matching
metadata = pd.read_csv("/home/jovyan/mount/gdrive/Pan_fetal/annotations/manifest_clean_120121.csv", index_col=0)
## Rename columns as they are in obs
metadata['donor'] = metadata['SAMPLE.NAME']
matching_df = metadata[["Sample.lanes","donor"]]
matching_df.columns = ["file","donor"]
matching_df[matching_df.file.isin(os.listdir(outdir))].to_csv("/nfs/team205/ed6/data/Fetal_immune/panfetal_file2donor.csv", header=False)