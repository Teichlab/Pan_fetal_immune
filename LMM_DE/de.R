library(Matrix)

bin_dir = "/nfs/team205/ed6/bin/Pan_fetal_immune/LMM_DE/"
input_data_dir = "/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_input_MYELOID/"
output_data_dir = "/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_MYELOID/"

source(paste0(bin_dir, "LMM.R"))

id=as.numeric(commandArgs()[5])

load(paste0(output_data_dir, "res.Rbin"))
load(paste0(output_data_dir, "Y_0.05_each_celltype.Rbin"))

ucelltypes=unique(matrix(unlist(strsplit(colnames(res$Z)[res$H[,names(res$nh)=="celltype_organ"]==1],":")),2)[1,])
for(i in ucelltypes[id]){
	de = getBFInt(Y,res,"celltype_organ",Celltype=i)
	saveRDS(de, file=paste(output_data_dir, "DE/de_",id,".RDS",sep=""))
}

