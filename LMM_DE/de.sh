#!/bin/bash

MEM=50000
if [ -z $1 ]
then
	bsub -J"atacJob[1-27]" -m modern_hardware -e /nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_MYELOID_PBULK/FarmOut/stderr.%I -o /nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_MYELOID_PBULK/FarmOut/stdout.%I -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team283 "sh de.sh \$LSB_JOBINDEX"
	exit
fi

/software/team170/miniconda3/bin/R --vanilla --quiet --args $1 < de.R


