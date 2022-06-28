#! /bin/bash
#BSUB -G cellgeni
#BSUB -J c2l.pred[1]
#BSUB -o %J.%I.c2l.pred.log
#BSUB -e %J.%I.c2l.pred.err
#BSUB -n 2
#BSUB -M32000
#BSUB -R "span[hosts=1] select[mem>32000] rusage[mem=32000]"
#BSUB -q gpu-normal
#BSUB -gpu "mode=shared:j_exclusive=yes:gmem=32000:num=1"

source activate cell2loc_env
cd /lustre/scratch117/cellgen/cellgeni/pasham/data/2202.c2l.service/tic-1349/prediction/alpha20
 
# v1 (tic-1349)
# vis=(/nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.spleen.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.spleen.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.spleen.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.thymus.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.thymus.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.thymus.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.liver.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.liver.h5ad 
#      /nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.liver.h5ad)

# redo one sample for TIC-1451
#the reference was taken from here /warehouse/cellgeni/tic-1349/ref/subsetSP/rsignatures/inf_aver.csv 
# vis=(/nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.batch1.spleen.h5ad)
# 
# outputs=(subsetSP 
#          subsetSP.maxAge13 
#          subsetSP.minAge14 
#          subsetTH 
#          subsetTH.maxAge13 
#          subsetTH.minAge14 
#          subsetLI 
#          subsetLI.maxAge13 
#          subsetLI.minAge14)
#redo one sample for TIC-1468
vis=(/nfs/team205/ed6/data/Fetal_immune/c2l_visium/fetal_immune_spatial.batch1.liver.h5ad)
outputs=(subsetLI)

/nfs/users/nfs_p/pm19/nfs/projects/2202.c2l.service/src/py/02.predict.cell.abundancies.py \
        ${vis[$LSB_JOBINDEX-1]} \
        ../../ref/${outputs[$LSB_JOBINDEX-1]}/rsignatures/inf_aver.csv \
        ${outputs[$LSB_JOBINDEX-1]}