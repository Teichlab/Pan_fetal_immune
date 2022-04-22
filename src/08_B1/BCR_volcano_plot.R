library(dplyr)
library(ggplot2)
library(ggrepel)
setwd('/home/jovyan/panfetal/')

results = read.csv('csv/b1_bcr_lr.csv',row.names = 1)
# add a column of NAs
results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
results$diffexpressed[results$logOddsRatio > 0.69 & results$adj_pval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
results$diffexpressed[results$logOddsRatio < -0.69 & results$adj_pval < 0.05] <- "DOWN"

results$delabel <- NA
results$delabel[results$diffexpressed != "NO"] <- rownames(results)[results$diffexpressed != "NO"]

pdf(paste0('/home/jovyan/mount/gdrive/Pan_fetal/plots_output/chenqu_jhub/',"volcano_plot.pdf"),width=6, height=4)
ggplot(data=results, aes(x=logOddsRatio, y=-log10(adj_pval), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#1F77B4", "black", "#D62728")) +
  xlab('log(odds ratio)')+
  xlim(-2,2)
dev.off()
