library(Matrix)

bin_dir = "/nfs/team205/ed6/bin/Pan_fetal_immune/LMM_DE/"
input_data_dir = "/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_input_MYELOID_PBULK/"
output_data_dir = "/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_MYELOID_PBULK/"

if (!dir.exists(output_data_dir)){
    system(paste("mkdir", output_data_dir))
}

source(paste0(bin_dir, "LMM.R"))
# count 
Y=t(readMM(paste0(input_data_dir,"matrix.mtx.gz"))) # original Y in mtx format
#load("Y.Rbin") # R binary of Y above (transposed)

# meta data
mdata=read.csv(paste0(input_data_dir,"metadata.csv.gz"),as.is=T)
names(mdata)[5]="celltype"

# celltype x organ interaction
celltype_organ=(paste(mdata$celltype,mdata$organ,sep=":"))
mdata=cbind(mdata,celltype_organ)

# cell proportion for each celltype
Z=model.matrix(~0+mdata$celltype)
P=(Y>0)%*%Z
P=t(t(P)/colSums(Z))
colnames(P)=gsub("mdata\\$celltype","",colnames(P))
save(P,file=paste0(output_data_dir, "P.Rbin"))

# gene used (prop of cell expressing the gene > 5%)
flag=apply(P>0.05,1,sum)>0
Y = Y[flag,]
save(flag,file=paste0(output_data_dir,"flag_0.05.Rbin"))
save(Y,P,Z,file=paste0(output_data_dir,"Y_0.05_each_celltype.Rbin"))

# N genes and N UMIs
mdata = cbind(mdata, n_gene=colSums(Y>0), n_umi=colSums(Y))

# linear mixed model
res=LFLMM(Y, mdata[,-c(1,2)], ITRMAX=300)
png(width=1000,height=1000,res=150,file=paste0(output_data_dir,"bar.png"));par(mar=c(7,3,1,1));Barplot(res);dev.off()

save(res,file=paste0(output_data_dir,"res.Rbin"))

# Prep for DE on celltypes
cmd <- paste0("mkdir ", output_data_dir, "DE")
system(cmd)
ucelltypes=unique(matrix(unlist(strsplit(colnames(res$Z)[res$H[,names(res$nh)=="celltype_organ"]==1],":")),2)[1,])
write.table(data.frame(ucelltypes), file=paste0(output_data_dir, "DE/celltype.txt"), col.names = FALSE)

# DE for marginal effects
de_celltype=getBF(Y,res,"celltype",DE1=1,AllClus=F)
save(de_celltype,file=paste0(output_data_dir, "DE/de_celltype.Rbin"))
de_organ=getBF(Y,res,"organ",DE1=1,AllClus=F)
save(de_organ,file=paste0(output_data_dir,"DE/de_organ.Rbin"))

