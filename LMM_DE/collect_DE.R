library(Matrix)

id = 'MYELOID_PBULK'
input_data_dir = paste0("/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_input_", id, "/")
output_data_dir = paste0("/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_output_", id, "/")

load(paste0(output_data_dir, "P.Rbin"))
load(paste0(output_data_dir,"flag_0.05.Rbin"))
celltype_table=read.table(paste0(output_data_dir,"DE/celltype.txt"),as.is=T)
gene=read.csv(paste0(input_data_dir, "gene.csv.gz"),as.is=T)[flag,]
mdata=read.csv(paste0(input_data_dir,"metadata.csv.gz"),as.is=T)
uorgan=names(table(mdata$organ))
# beta = ltsr = delta = gname = nct = NULL
th = 0.9
res_df = data.frame()
for(celltype in celltype_table[,1]){
	# DE result
	de=readRDS(paste(output_data_dir, "DE/de_",celltype,".RDS",sep=""))
    # Exclude if celltype is present just in one organ
    if(dim(de$beta)[2] > 1){
        # gene is expressed by at least 5% of cells
        flag5 = P[flag,celltype]>0.05

        # DE genes i.e., LTSR>th in any contrast
        flag_de = apply(de$ltsr>th,1,sum)>0 & flag5

        # average expression
        rownames(de$beta) = 1:nrow(de$beta) # replaced to be gene names
        beta1 = de$beta[flag_de,match(uorgan, matrix(unlist(strsplit(colnames(de$beta),":")),2)[2,])]

        # ltsr
        ltsr1 = de$ltsr[flag_de,match(uorgan, gsub(":","",substring(colnames(de$ltsr),nchar(colnames(de$ltsr))-2,nchar(colnames(de$ltsr)))))]

        # log fold change
        delta1 = de$deltabeta[flag_de,match(uorgan, gsub(":","",substring(colnames(de$ltsr),nchar(colnames(de$ltsr))-2,nchar(colnames(de$ltsr)))))]

        # ordering
        imax = function(x){seq(x)[x%in%max(x,na.rm=T)][1]}
        maxorg = apply(ltsr1,1,imax)
        D1 = diag(length(uorgan))[maxorg,] # masking 
        d1norm = (rank(delta1)/length(c(delta1)))*D1
        d1norm[is.na(d1norm)]=0
        ord1 = order(-(colSums(t(d1norm) + t(D1)*rev(seq(length(uorgan))))))

        # collect
        colnames(beta1) = colnames(delta1) = colnames(ltsr1) = uorgan
        gname = gene[as.numeric(rownames(beta1)),2][ord1]
        beta = beta1[ord1,]
        ltsr = ltsr1[ord1,]
        delta = delta1[ord1,]
        nct = nrow(beta1)
    #         beta = rbind(beta, beta1[ord1,])
    #         ltsr = rbind(ltsr, ltsr1[ord1,])
    #         delta = rbind(delta, delta1[ord1,])
    #         nct = c(nct, nrow(beta1))   
        rownames(beta)=gname
        morgan = uorgan[apply(ltsr,1,imax)]
        mdelta = apply(delta*(ltsr>th),1,sum)
#         print(head(delta))
#         print(head(morgan))
        mltsr = apply(ltsr,1,max,na.rm=T)
        df = data.frame(celltype=rep(celltype_table[celltype,2],nct),organ=morgan,gene=gname,ltsr=mltsr,logF=mdelta,beta)
        res_df = rbind(res_df,df)
    }
}

## Save
write.csv(res_df, file=paste0(output_data_dir,"DE/DE_summary_", id, "_lstr", th, ".csv"))
