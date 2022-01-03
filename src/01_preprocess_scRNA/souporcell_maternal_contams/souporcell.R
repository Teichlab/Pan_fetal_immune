setwd('~/panfetal/souporcell/')
library(ggplot2)

annot = read.csv('/nfs/team205/ed6/data/Fetal_immune/scVI_outs/PAN.A01.v01.entire_data_raw_count.20210429.scVI_out.clustering.csv')
#annot[startsWith(annot$file,'Human'),'index'] = paste0(annot[startsWith(annot$file,'Human'),'index'],'-1')

## build the data frame
df = data.frame(row.names = c('k=1','k=2','k=3','log_prob improvement','BIC_k=1','BIC_k=2','BIC_k=3','perc_unassigned_doublet','% smaller genotype','leiden_150_no_small_genotype','leiden_150_no_big_genotype'))
BIC = data.frame(row.names = c('loci_1','loci_2','loci_3','no_cells'), stringsAsFactors = FALSE)
path='/lustre/scratch117/cellgen/team205/sharedData/cs42/panfetal-donor-souporcell/'
donors = list.dirs(path = path, full.names = FALSE, recursive = FALSE)

## read log_probabilities when k=1,2,3, and calculate (k=3-k=2)/(k=2-k=1)
for (j in 1:length(donors)){
  donor = donors[j]
  for (i in 1:3){
    k=i
    path_new = paste0(path, donor, '/',k,'/clusters.err')
    err = readLines(path_new)
    log_prob = as.numeric(strsplit(err[length(err)], split='=')[[1]][2])
    df[k,donor]= log_prob
    
    loci_used = as.numeric(strsplit(err[1], split=' ')[[1]][4])
    BIC[k,donor] = loci_used
  }

  df[4, donor] = (df[3, donor] - df[2, donor]) / (df[2, donor] - df[1, donor]) # calculate (k=3-k=2)/(k=2-k=1)
  
  ## read in souporcell output
  souporcell_output = read.csv(paste0(path, donor, '/','2','/clusters.tsv'), sep='\t')
  
  ## calculate BIC scores
  # no. of cells in that donor
  BIC[4, donor] = as.numeric(nrow(souporcell_output))
  for (i in 1:3){
    k = i 
    df[i+4, donor] = -2*df[i,donor] + i*BIC[k,donor]*log(BIC[4, donor])
  }
  
  ## calculate % (unassigned + doublet)
  df['perc_unassigned_doublet', donor] = sum(souporcell_output$status %in% c('unassigned','doublet')) / nrow(souporcell_output)
  
  ## calculate % of smaller genotype
  df['% smaller genotype', donor] = sum(souporcell_output$assignment =='0') / (sum(souporcell_output$assignment =='0') + sum(souporcell_output$assignment =='1'))
  if (df['% smaller genotype', donor] >0.5){
    df['% smaller genotype', donor] = 1 - df['% smaller genotype', donor]
  }
  
  ## number of unique leiden_150 in the smaller genotype cluster
  # unify the barcodes first
  for (l in 1:nrow(souporcell_output)){
    barcode = souporcell_output$barcode[l]
    split = strsplit(barcode,split='-')[[1]]
    barcode_new = paste0(split[1],'-',split[2])
    souporcell_output$barcode[l] = barcode_new
  }
  
  # lookup the leiden_150
  souporcell_output$annot = annot$leiden_150[match(souporcell_output$barcode,annot$X)]
  souporcell_output_0 = souporcell_output[souporcell_output$assignment == '0',]
  no0 = length(unique(souporcell_output_0$annot[!(souporcell_output_0$annot %in% c('',NA))]))
  souporcell_output_1 = souporcell_output[souporcell_output$assignment == '1',]
  no1 = length(unique(souporcell_output_1$annot[!(souporcell_output_1$annot %in% c('',NA))]))
  
  if (sum(souporcell_output$assignment =='0') < sum(souporcell_output$assignment =='1')){
    df['leiden_150_no_small_genotype',donor] = no0
    df['leiden_150_no_big_genotype',donor] = no1
  }
  if (sum(souporcell_output$assignment =='0') > sum(souporcell_output$assignment =='1')){
    df['leiden_150_no_small_genotype',donor] = no1
    df['leiden_150_no_big_genotype',donor] = no0
  }
}

# Is BIC_k=2 the smallest
n_donor = ncol(df)

df['BIC_k=2 smallest',] = FALSE
for (i in 1:n_donor){
  df['BIC_k=2 smallest',i] = df['BIC_k=2',i] < df['BIC_k=1',i] && df['BIC_k=2',i] < df['BIC_k=3',i]
}

df['BIC_k=3 smallest',] = FALSE
for (i in 1:n_donor){
  df['BIC_k=3 smallest',i] = df['BIC_k=3',i] < df['BIC_k=2',i] && df['BIC_k=3',i] < df['BIC_k=1',i]
}

# plot BIC scores
df_selected = df[,df['BIC_k=2 smallest',]==1] # F19, F37
n_donor_selected = ncol(df_selected)
df2 = data.frame('metric_name' = c(rep('k=1',n_donor_selected), rep('k=2',n_donor_selected),rep('k=3',n_donor_selected)),
                 'metric' = c(as.numeric(df_selected['BIC_k=1',]), as.numeric(df_selected['BIC_k=2',]),as.numeric(df_selected['BIC_k=3',])),
                 'donor' = rep(colnames(df_selected),3))
ggplot(data = df2, mapping = aes(y = metric, x = metric_name, group = donor, colour=donor)) +
  geom_line(size=1.2) +
  geom_point(size=5) +
  labs(y= "BIC score", x = "k")


# plot log_prob improvement & %(unassigned+doublet)
df_new = data.frame('metric_name' = c(rep('log_prob improvement',n_donor), rep('%(unassigned+doublet)',n_donor)),
                    'metric' = c(as.numeric(df[4,]), as.numeric(df[5,])),
                    'donor' = rep(c(1:n_donor),2))

ggplot(df_new[c(1:25),], aes(x=metric_name,y=metric))+
  geom_point(color = as.numeric(df_new$donor))+
  scale_color_brewer(palette="Set3")


df3 = data.frame('k'=rownames(df)[1:3], 'log_prob'=df$F19[1:3])
ggplot(df3, aes(x=k,y=log_prob))+
  geom_point()

write.csv(t(df), file = 'souporcell_output.csv')

