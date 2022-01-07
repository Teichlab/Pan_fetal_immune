#!/bin/bash
set -e pipefail

PATH="/software/singularity-v3.6.4/bin:$PATH"

META=/lustre/scratch117/cellgen/team205/cs42/panfetal/entire_data_cpdb_meta.tsv
COUNTS=/lustre/scratch117/cellgen/team205/cs42/panfetal/entire_data_cpdb.h5ad
OUTPUT_PATH=/lustre/scratch117/cellgen/team205/sharedData/cs42/panfetal-cellphonedb/output

singularity run --bind /nfs,/lustre  \
  /nfs/cellgeni/singularity/images/cellphonedb-v3.0.0-alpha.sif \
  method statistical_analysis \
  $META $COUNTS \
  --output-path=$OUTPUT_PATH \
  --counts-data=hgnc_symbol