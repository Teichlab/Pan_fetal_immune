#!/bin/bash

source activate emma_env
cd /nfs/team205/ed6/bin/Pan_fetal_immune/

## Split 1
python ./utils/make_splits.py PFI

## Split 2
python ./utils/cluster_split.py STROMA
python ./utils/cluster_split.py MEM_PROGENITORS
python ./utils/cluster_split.py MYELOID_LYMPHOID

python ./utils/make_splits.py MYELOID_LYMPHOID

## Split 3
python ./utils/cluster_split.py MYELOID
python ./utils/cluster_split.py LYMPHOID

python ./utils/make_splits.py LYMPHOID

## Split 4
python ./utils/cluster_split.py NKT
python ./utils/cluster_split.py B