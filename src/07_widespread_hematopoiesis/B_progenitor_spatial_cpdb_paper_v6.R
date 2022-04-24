setwd('/home/jovyan/panfetal')
library("dplyr")   
library(reshape2)

B_prog = c('PRE_PRO_B','PRO_B','LATE_PRO_B','LARGE_PRE_B','SMALL_PRE_B')
organs = c('LI','TH','SP')

#cpdb_path = '/lustre/scratch117/cellgen/team205/sharedData/cs42/cellphonedb/organs/'
cpdb_path = '/lustre/scratch117/cellgen/team205/sharedData/cs42/cellphonedb/organs2/' # with updated visium mapping

#### from Fig 3D, co-localising cell types across liver/spleen/thymus
other_players = c('NK', 'CYCLING_NK', 'MACROPHAGE_LYVE1_HIGH', 'TYPE_1_INNATE_T', 'ILC3') # didn't include LMPP_MLP as likely upstream progenitor

########################## look in cpdb #############################
cpdb=NULL

for (organ in organs){
  means = read.csv(paste0(cpdb_path, organ, '/means.txt'), sep = '\t',stringsAsFactors = F)
  pvals = read.csv(paste0(cpdb_path, organ, '/pvalues.txt'), sep = '\t',stringsAsFactors = F)
  
  ligand_pairs <- t(sapply(strsplit(means$interacting_pair, '_'), function(v){c(v[1], paste(v[-1], collapse='.'))}))
  ct_pairs <- t(sapply(strsplit(colnames(means), '[.]'), function(v){c(v[1], paste(v[-1], collapse='.'))}))
  
  column_select <- paste(ct_pairs[,1], ct_pairs[,2]) %in% c(outer(B_prog, other_players, paste), outer(other_players, B_prog, paste))
  means_selected <- means[, column_select]
  pvals_selected <- pvals[, column_select]
  interacting_pair <- paste0(ligand_pairs[,1], '|', ligand_pairs[,2])
  colnames(means_selected) <- colnames(pvals_selected) <- 
    paste0(ct_pairs[column_select,1], '|', ct_pairs[column_select,2])
  means_selected <- cbind(interacting_pair, means_selected)
  pvals_selected <- cbind(interacting_pair, pvals_selected)
  
  # aggregate means
  means_selected_sub <- means_selected[apply(pvals_selected[, -1], 1, min) < 0.05, ]
  df <- melt(means_selected_sub, id='interacting_pair')
  df$variable = as.character(df$variable)
  
  df$cell_1 <- sapply(strsplit(df$variable, '|', fixed='TRUE'), function(v){v[1]})
  df$cell_2 <- sapply(strsplit(df$variable, '|', fixed='TRUE'), function(v){v[2]})
  df$ligand_1 <- sapply(strsplit(df$interacting_pair, '|', fixed='TRUE'), function(v){v[1]})
  df$ligand_2 <- sapply(strsplit(df$interacting_pair, '|', fixed='TRUE'), function(v){v[2]})
  to_swap <- !(df$cell_1 %in% B_prog) # put B_progenitor as the first cell
  df[to_swap, c(4,5,6,7)] <- df[to_swap, c(5,4,7,6)]
  df$cell_pair <- paste0(df$cell_1, '.', df$cell_2)
  df$ligand_pair <- paste0(df$ligand_1, '_', df$ligand_2)
  df_new = aggregate(value~cell_2+ligand_pair, data=df, FUN=mean) # aggregating all B progenitors for the same interacting cell & ligand_pair
  means_update = dcast(df_new, ligand_pair ~ cell_2)
  
  # aggregate pvals
  pvals_selected_sub <- pvals_selected[apply(pvals_selected[, -1], 1, min) < 0.05, ]
  df <- melt(pvals_selected_sub, id='interacting_pair')
  df$variable = as.character(df$variable)
  
  df$cell_1 <- sapply(strsplit(df$variable, '|', fixed='TRUE'), function(v){v[1]})
  df$cell_2 <- sapply(strsplit(df$variable, '|', fixed='TRUE'), function(v){v[2]})
  df$ligand_1 <- sapply(strsplit(df$interacting_pair, '|', fixed='TRUE'), function(v){v[1]})
  df$ligand_2 <- sapply(strsplit(df$interacting_pair, '|', fixed='TRUE'), function(v){v[2]})
  to_swap <- !(df$cell_1 %in% B_prog)
  df[to_swap, c(4,5,6,7)] <- df[to_swap, c(5,4,7,6)]
  df$cell_pair <- paste0(df$cell_1, '.', df$cell_2)
  df$ligand_pair <- paste0(df$ligand_1, '_', df$ligand_2)
  df_new = aggregate(value~cell_2+ligand_pair, data=df, FUN=min) # aggregating all B progenitors for the same interacting cell & ligand_pair
  pvals_update = dcast(df_new, ligand_pair ~ cell_2)
  
  colnames(means_update)[-1] <- paste0(organ, '-', colnames(means_update)[-1])
  colnames(pvals_update)[-1] <- paste0(organ, '-', colnames(pvals_update)[-1])
  cpdb[[organ]]=list(means_update, pvals_update)
}


# find intersection of ligand_pair across organs
merged_means <- cpdb[['LI']][[1]]
merged_pvals <- cpdb[['LI']][[2]]

for (organ in c('SP','TH')){
  merged_means <- merge(merged_means, cpdb[[organ]][[1]], by='ligand_pair', all.x=FALSE)
  merged_pvals <- merge(merged_pvals, cpdb[[organ]][[2]], by='ligand_pair', all.x=FALSE)
}

# only display ligand_pair that are significant in three organs & rank by maximum mean
scores <- apply(merged_means[,-1], 1, max) * (apply(merged_pvals[, -1], 1, 
                                                    function(v){max(min(v[c(1,2,3,4,5)]), min(v[c(6,7,8,9,10)]), min(v[c(11,12,13,14,15)]))}) < 0.05)

merged_means_filtered = head(merged_means[order(scores, decreasing=TRUE), ], 60)
merged_pvals_filtered = head(merged_pvals[order(scores, decreasing=TRUE), ], 60)

write_csv(merged_means_filtered, file='csv/B_prog_cpdb_means_filtered.csv')
write_csv(merged_pvals_filtered, file='csv/B_prog_cpdb_pvals_filtered.csv')
