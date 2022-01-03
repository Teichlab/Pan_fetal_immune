setwd('~/panfetal/souporcell/')

path='/lustre/scratch117/cellgen/team205/sharedData/cs42/panfetal-donor-souporcell/'
donors = c('F19','F33','F37')
maternal_barcodes = c()

for (i in 1:3){
  donor = donors[i]
  
  ## read in souporcell output
  souporcell_output = read.csv(paste0(path, donor, '/','2','/clusters.tsv'), sep='\t')
  
  if (sum(souporcell_output$assignment==0) > sum(souporcell_output$assignment==1)){
    maternal_genotype = 1
  }
  else {
    maternal_genotype = 0
  }
  maternal_barcodes = c(maternal_barcodes, souporcell_output$barcode[souporcell_output$assignment==maternal_genotype])
}

write.csv(maternal_barcodes, file = 'maternal_barcodes.csv')

