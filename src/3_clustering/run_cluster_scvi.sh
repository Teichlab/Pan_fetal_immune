#!/bin/bash

MEM=25000
input_dir=/nfs/team205/ed6/data/Fetal_immune/scVI_outs/

emb_files=$(ls $input_dir/*ldims.npy)

for f in $emb_files
    do 
    out_prefix=$(echo $f | cut -d/ -f 9 | sed 's/PAN.A01.v01.entire_data_raw_count.wGut.//' | sed 's/.npy//')
    echo "python PFI_splits_clustering/cluster_scvi.py $f" | \
    bsub -e $input_dir/${out_prefix}.err -o $input_dir/${out_prefix}.out -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team283
    done
