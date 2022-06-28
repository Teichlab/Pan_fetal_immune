#! /bin/bash
#BSUB -G cellgeni
#BSUB -J c2l.ref[1-8]
#BSUB -o %J.%I.c2l.ref.log
#BSUB -e %J.%I.c2l.ref.err
#BSUB -n 2
#BSUB -M32000
#BSUB -R "span[hosts=1] select[mem>32000] rusage[mem=32000]"
#BSUB -q gpu-normal
#BSUB -gpu "mode=shared:j_exclusive=yes:gmem=32000:num=1"

source activate cell2loc_env
cd /nfs/users/nfs_p/pm19/nfs/projects/2202.c2l.service/out/ref


inputs=(/nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetSP.exclude_lowQ.h5ad
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetSP.maxAge13.exclude_lowQ.h5ad
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetSP.minAge14.exclude_lowQ.h5ad 
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetTH.exclude_lowQ.addTECs.keep_fetal_TECs.h5ad 
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetTH.maxAge13.exclude_lowQ.addTECs.keep_fetal_TECs.h5ad 
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetTH.minAge14.exclude_lowQ.addTECs.keep_fetal_TECs.h5ad 
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetLI.exclude_lowQ.h5ad 
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetLI.maxAge13.exclude_lowQ.h5ad 
        /nfs/team205/ed6/data/Fetal_immune/c2l_scRNA_references/PAN.A01.v01.c2l_reference.v2.subsetLI.minAge14.exclude_lowQ.h5ad)

outputs=(subsetSP 
         subsetSP.maxAge13 
         subsetSP.minAge14 
         subsetTH 
         subsetTH.maxAge13 
         subsetTH.minAge14 
         subsetLI 
         subsetLI.maxAge13 
         subsetLI.minAge14)

../../src/py/01.estimate.signatures.py ${inputs[$LSB_JOBINDEX-1]} ${outputs[$LSB_JOBINDEX-1]}