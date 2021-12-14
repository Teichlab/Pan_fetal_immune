setwd('/home/jovyan/panfetal')
library("dplyr")   
library(reshape2)

B_prog = c('PRE_PRO_B','PRO_B','LATE_PRO_B','LARGE_PRE_B','SMALL_PRE_B')
organs = c('LI','TH','SP')

#cpdb_path = '/lustre/scratch117/cellgen/team205/sharedData/cs42/cellphonedb/organs/'
cpdb_path = '/lustre/scratch117/cellgen/team205/sharedData/cs42/cellphonedb/organs2/' # with updated visium mapping

#### from Fig 3D, 3 co-localising cell types across liver/spleen/thymus
other_players = c('ILC3','MACROPHAGE_LYVE1_HIGH','NK')

########################## look in cpdb #############################
cpdb=NULL

for (organ in organs){
  sig_means = read.csv(paste0(cpdb_path, organ, '/significant_means.txt'), sep = '\t',stringsAsFactors = F)
  
  cpdb[[organ]]= sig_means
  
  list_select = rep(NA, ncol(sig_means)-12)
  for (i in 1:ncol(sig_means)-12){
    cell_A = strsplit(colnames(sig_means)[i+12], split = '[.]')[[1]][1]
    cell_B = strsplit(colnames(sig_means)[i+12], split = '[.]')[[1]][2]
    
    if (cell_A %in% B_prog){
      list_select[i] = cell_B %in% other_players
    } else if (cell_A %in% other_players){
      list_select[i] = cell_B %in% B_prog
    } else {
      list_select[i] = FALSE
    }
  }
  
  sig_selected = sig_means[,c(FALSE, TRUE, rep(FALSE, 10),list_select)] #only keep 'interacting_pair' and the selected cell type pairs
  
  # change all NA to 0
  sig_selected[is.na(sig_selected)] <- 0
  
  # filter for rows with mean expression >0  
  row_sum = rep(0,nrow(sig_selected))
  for (i in 1:nrow(sig_selected)){
    row_sum[i] = sum(sig_selected[i,c(2:ncol(sig_selected))] > 0) 
  }
  sig_selected_sub = sig_selected[row_sum >0,]
  
  #rownames(sig_selected_sub) = sig_selected_sub$interacting_pair
  #sig_selected_sub = sig_selected_sub[,2:ncol(sig_selected_sub)]
  
  ### change order of cell types to B_prog first
  sig_selected_sub = sig_selected[row_sum >0,]
  df = melt(sig_selected_sub, id=c('interacting_pair'))
  df$interacting_pair = as.character(df$interacting_pair)
  df$variable = as.character(df$variable)
  df$cell_1 <- sapply(strsplit(df$variable, '.', fixed='TRUE'), function(v){v[1]})
  df$cell_2 <- sapply(strsplit(df$variable, '.', fixed='TRUE'), function(v){v[2]})
  df$ligand_1 <- sapply(strsplit(df$interacting_pair, '_'), function(v){v[1]})
  df$ligand_2 <- sapply(strsplit(df$interacting_pair, '_'), function(v){v[2]})
  to_swap <- !(df$cell_1 %in% B_prog)
  df[to_swap, c(4,5,6,7)] <- df[to_swap, c(5,4,7,6)]
  df$cell_pair <- paste0(df$cell_1, '.', df$cell_2)
  df$ligand_pair <- paste0(df$ligand_1, '_', df$ligand_2)
  df_new = aggregate(value~cell_2+ligand_pair, data=df, FUN=mean) # aggregating all B progenitors for the same interacting cell & ligand_pair
  df_update = dcast(df_new, ligand_pair ~ cell_2)
  #rownames(df_new) = df_new$ligand_pair
  #df_new=df_new[,2:3]
  df_update = df_update[rowSums(df_update[,2:ncol(df_update)])>0,]
  
  colnames(df_update)[-1] = paste0(organ,'-',colnames(df_update))[-1]
  cpdb[[organ]]=df_update
}

# find intersection of ligand_pair across organs
merged_list <- cpdb[['LI']]
for (organ in c('SP','TH')){
  merged_list <- merge(merged_list, cpdb[[organ]], by='ligand_pair', all.x=FALSE)
}

merged_list_filtered = merged_list[apply(merged_list[,-1], 1, max) > 1,]

write_csv(merged_list_filtered, file='csv/B_prog_cpdb_all.csv')
