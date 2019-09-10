## 04c - Heatmap for differentially expressed genes between control and injury
# Plotting done with the results of differential expression in R between SST-RMA normalized signals using all 16 samples (both FC and LA; all three conditions together).
# Only LA samples are plotted here

setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/injury_DE_res.rda', verbose=T)
library(jaffelab)
library(ggplot2)
library(ggrepel)
library(limma)
theme_set(theme_bw(base_size=40) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))

## subsetting to plot only FDR significant CTRL v INJURY trasncript clusters in LA and LA samples

signal_subset = sst_rma_signal[as.character(merged_SstRma_disease[merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-adj.P.Val`<0.10,'transcriptCluster']),as.character(pd[pd$brain_region=="LA",'sampleNames'])]
  
##
library(pheatmap)
library(RColorBrewer)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
##

annotation <- data.frame( Dx = pd[pd$brain_region=="LA",'Group'] )

rownames(annotation) <- pd[pd$brain_region=="LA",'sampleNames']

#ann_colors = list(Dx = c("black","#ab1323" ) )
#names(ann_colors[['Dx']])<- levels(annotation$Dx)

  # list with colors for each annotation
  mat_colors <- list(group=brewer.pal(3,"Set1"))
  names(mat_colors$group) = unique(pd$Group)
  names(mat_colors)<-"Dx"
  
library(viridis)
pdf('plots/heatmap_amygdalaOnly_FDR10_CTRL_v_InjuryGenes.pdf',height=12,width=12,onefile=F,useDingbats=FALSE)
pheatmap::pheatmap(signal_subset,
                   cluster_cols = FALSE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   treeheight_row=100,
                   treeheight_col=50,
                   cutree_cols = 2,
                   gaps_col = c(3,5),
                   labels_row = as.character(merged_SstRma_disease[rownames(signal_subset),'Symbol']),
                   annotation_col = annotation,
                   annotation_colors=mat_colors,
                   show_colnames = FALSE, 
                   show_rownames = TRUE,
                   fontsize = 10,
                   fontsize_row=5,
                   scale='row',
                   width=12,
                   height=12,
                   color = inferno(20),
                   legend_labels = c('Control', 'Injury', 'Injury + Tx') )
dev.off()


tiff('plots/heatmap_amygdalaOnly_FDR10_CTRL_v_InjuryGenes.tiff', height=12,width=12, units = "in", res = 600)
pheatmap::pheatmap(signal_subset,
                   cluster_cols = FALSE,
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   treeheight_row=100,
                   treeheight_col=50,
                   cutree_cols = 2,
                   gaps_col = c(3,5),
                   labels_row = as.character(merged_SstRma_disease[rownames(signal_subset),'Symbol']),
                   annotation_col = annotation,
                   annotation_colors=mat_colors,
                   show_colnames = FALSE, 
                   show_rownames = TRUE,
                   fontsize = 10,
                   fontsize_row=5,
                   scale='row',
                   width=12,
                   height=12,
                   color = inferno(20),
                   legend_labels = c('Control', 'Injury', 'Injury + Tx') )
dev.off()


