#!/bin/bash

MEM=50000
if [ -z $1 ]
then
	bsub -q long -m modern_hardware -e /nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_MYELOID_PBULK/stderr -o /nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_MYELOID_PBULK/stdout -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team283 "sh run.sh farm"
	exit
fi

/software/team170/miniconda3/envs/seurat/bin/R --vanilla --quiet --args < run.R


