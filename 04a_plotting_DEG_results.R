## 04a - Plotting differential expression effects between groups.
# Plotting done with the results of differential expression in R between SST-RMA normalized signals using all 16 samples (both FC and LA; all three conditions together).

setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/normalized_data_n16_sst-rma.rda', verbose=T)
load('rda/FullModel_mergedStats_sst-RMA_R.rda', verbose=T)
library(jaffelab)
library(ggplot2)
library(ggrepel)
library(limma)
theme_set(theme_bw(base_size=40) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))

merged_SstRma_disease = mergedStats
#############====== Control vs. PYR+LPS =======################
merged_SstRma_disease = dplyr::mutate(as.data.frame(merged_SstRma_disease), LA_sig=ifelse(merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-adj.P.Val`< 0.10, TRUE, FALSE))
merged_SstRma_disease = dplyr::mutate(as.data.frame(merged_SstRma_disease), FC_sig=ifelse(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-adj.P.Val`< 0.10, TRUE, FALSE))
merged_SstRma_disease$Sig = "Not sig"
merged_SstRma_disease$Sig[merged_SstRma_disease$FC_sig] = "FC Sig"
merged_SstRma_disease$Sig[merged_SstRma_disease$LA_sig] = "LA Sig"
merged_SstRma_disease$Sig[merged_SstRma_disease$LA_sig &  merged_SstRma_disease$FC_sig] = "Both Sig"
merged_SstRma_disease$Sig = as.factor(merged_SstRma_disease$Sig)

table(merged_SstRma_disease$Sig )

rownames(merged_SstRma_disease) = merged_SstRma_disease$transcriptCluster

##### Scatterplot of t-statistics between regions
lfc_comp = ggplot(data=merged_SstRma_disease, aes(x=`FC_CTRL_v_PYR_LPS-t`, y= `LA_CTRL_v_PYR_LPS-t`))				 
lfc_comp = lfc_comp + 
           geom_point() + 
           labs(x="Frontal cortex T-statistic", y="Lateral amygdala T-statistic", title = 'CTRL (N=3) vs. PYR+LPS (N=2)') + 
           geom_hline(yintercept=0,colour='red',size=1) + 
           geom_vline(xintercept=0,colour='red',size=1)
          # theme(legend.position='none') +
          # scale_colour_manual(values=c("red","steelblue", "green3","grey") ) + 
#genes_of_high_interest = merged_SstRma_disease[(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC`>1.5 & merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`>1.5) | (merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC` < -1.5 & merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`< -1.5),]
#lfc_comp = lfc_comp + geom_text_repel(data=genes_of_high_interest, aes(label=`Symbol`))

ggsave(lfc_comp, file = "plots/Control vs PYR+LPS T-statistics.pdf",
       height=12,width=12,units="in")		   

ggsave(lfc_comp, file = "plots/Control vs PYR+LPS T-statistics.tiff",
       height=12,width=12,units="in", dpi=600)		  

cor.test(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-t`, merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-t`)
#cor.test(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC`, merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`)

##### Add in average expression values (for condition-region, signal post-sstRma normalization w/ log transform)
merged_SstRma_disease$`meanSignal_FC-Control` <- rowMeans(sst_rma_signal[rownames(merged_SstRma_disease),pd$sampleNames[pd$Group=="CTRL"&pd$brain_region=="FC"] ])
merged_SstRma_disease$`meanSignal_FC-Injury` <- rowMeans(sst_rma_signal[rownames(merged_SstRma_disease),pd$sampleNames[pd$Group=="PYR_LPS"&pd$brain_region=="FC"] ])

merged_SstRma_disease$`meanSignal_LA-Control` <- rowMeans(sst_rma_signal[rownames(merged_SstRma_disease),pd$sampleNames[pd$Group=="CTRL"&pd$brain_region=="LA"] ])
merged_SstRma_disease$`meanSignal_LA-Injury` <- rowMeans(sst_rma_signal[rownames(merged_SstRma_disease),pd$sampleNames[pd$Group=="PYR_LPS"&pd$brain_region=="LA"] ])


merged_SstRma_disease$Region_MinFDR = matrixStats::rowMins(as.matrix(merged_SstRma_disease[,c('FC_CTRL_v_PYR_LPS-adj.P.Val', 'LA_CTRL_v_PYR_LPS-adj.P.Val')]) )
merged_SstRma_disease$Region_MinP = matrixStats::rowMins(as.matrix(merged_SstRma_disease[,c('FC_CTRL_v_PYR_LPS-P.Value', 'LA_CTRL_v_PYR_LPS-P.Value')]) )


save(merged_SstRma_disease, sst_rma_signal, pd, file='rda/injury_DE_res.rda')
##### Table of FDR significant changes between control and injury model

merged_SstRma_disease_save = merged_SstRma_disease[merged_SstRma_disease$Sig!="Not sig", ]

merged_SstRma_disease_save$regionConsistent = FALSE
merged_SstRma_disease_save$regionConsistent[merged_SstRma_disease_save$Sig=="Both Sig"] <- TRUE

merged_SstRma_disease_save$regionConsistent[merged_SstRma_disease_save$Sig=="FC Sig" & merged_SstRma_disease_save$`LA_CTRL_v_PYR_LPS-P.Value`<0.05 & sign(merged_SstRma_disease_save[, 'FC_CTRL_v_PYR_LPS-logFC']) == sign(merged_SstRma_disease_save[, 'LA_CTRL_v_PYR_LPS-logFC'])] <- TRUE

merged_SstRma_disease_save$regionConsistent[merged_SstRma_disease_save$Sig=="LA Sig" & merged_SstRma_disease_save$`FC_CTRL_v_PYR_LPS-P.Value`<0.05 & sign(merged_SstRma_disease_save[, 'LA_CTRL_v_PYR_LPS-logFC']) == sign(merged_SstRma_disease_save[, 'FC_CTRL_v_PYR_LPS-logFC'])] <- TRUE

  
write.csv(merged_SstRma_disease_save[order(merged_SstRma_disease_save$Region_MinFDR, merged_SstRma_disease_save$Region_MinP),c('transcriptCluster','Symbol','EntrezID','Accession','gene_description','Sig','regionConsistent','Region_MinFDR','Region_MinP','FC_CTRL_v_PYR_LPS-adj.P.Val', 'LA_CTRL_v_PYR_LPS-adj.P.Val','FC_CTRL_v_PYR_LPS-P.Value', 'LA_CTRL_v_PYR_LPS-P.Value','FC_CTRL_v_PYR_LPS-logFC', 'LA_CTRL_v_PYR_LPS-logFC', 'meanSignal_FC-Control','meanSignal_FC-Injury','meanSignal_LA-Control', 'meanSignal_LA-Injury')],file='csvs/CTRL_v_PYR+LPS_FDR_10.csv',row.names=FALSE)


### Check overlap stats for genes that are implicated in one region but not supported at FDR in the other
#Frontal cortex genes
table( sign(merged_SstRma_disease[merged_SstRma_disease$Sig=="Both Sig", 'FC_CTRL_v_PYR_LPS-logFC']) == sign(merged_SstRma_disease[merged_SstRma_disease$Sig=="Both Sig", 'LA_CTRL_v_PYR_LPS-logFC'])  )


table(merged_SstRma_disease[merged_SstRma_disease$Sig=="FC Sig", 'LA_CTRL_v_PYR_LPS-P.Value']<0.05 & sign(merged_SstRma_disease[merged_SstRma_disease$Sig=="FC Sig", 'FC_CTRL_v_PYR_LPS-logFC']) == sign(merged_SstRma_disease[merged_SstRma_disease$Sig=="FC Sig", 'LA_CTRL_v_PYR_LPS-logFC'])  )

#Lateral amygdala genes
table(merged_SstRma_disease[merged_SstRma_disease$Sig=="LA Sig", 'FC_CTRL_v_PYR_LPS-P.Value']<0.05 & sign(merged_SstRma_disease[merged_SstRma_disease$Sig=="LA Sig", 'LA_CTRL_v_PYR_LPS-logFC']) == sign(merged_SstRma_disease[merged_SstRma_disease$Sig=="LA Sig", 'FC_CTRL_v_PYR_LPS-logFC'])  )


## get some more stats
length(merged_SstRma_disease_save$transcriptCluster)
length(unique(merged_SstRma_disease_save$Symbol))
table(merged_SstRma_disease_save$Sig)

#####  Volcano Facet Plot
LA_max_P_FDR10 = max(merged_SstRma_disease[merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-adj.P.Val`<.10,'LA_CTRL_v_PYR_LPS-P.Value'])
FC_max_P_FDR10 = max(merged_SstRma_disease[merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-adj.P.Val`<.10,'FC_CTRL_v_PYR_LPS-P.Value'])
shared_thresh = min(FC_max_P_FDR10,LA_max_P_FDR10)

LA_results = dplyr::mutate(as.data.frame(merged_SstRma_disease), sig=ifelse(merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-adj.P.Val`< 0.10, "FDR<0.10", "FDR>0.10"))
LA_results <- LA_results[,c("Symbol","LA_CTRL_v_PYR_LPS-P.Value","LA_CTRL_v_PYR_LPS-logFC","sig")]
colnames(LA_results) = c("GeneSymbol", "P-value", "log2FC", "Sig")
LA_results$Type = "Lateral Amygdala"

FC_results = dplyr::mutate(as.data.frame(merged_SstRma_disease), sig=ifelse(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-adj.P.Val`< 0.10, "FDR<0.10", "FDR>0.10"))
FC_results <- FC_results[,c("Symbol","FC_CTRL_v_PYR_LPS-P.Value","FC_CTRL_v_PYR_LPS-logFC","sig")]
colnames(FC_results) = c("GeneSymbol", "P-value", "log2FC", "Sig")
FC_results$Type = "Frontal Cortex"

all_results <- rbind(LA_results, FC_results)

## Plotting
all_volc <- ggplot(all_results, aes(log2FC, -log10(`P-value`))) +
  geom_point(aes(col=Sig)) +
  facet_wrap(~ Type, ncol= 2, nrow=1) +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"(pvalue)")) + 
  scale_color_manual(values=c("red", "grey")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, round(max(-log10(all_results$`P-value`))+.5) ) ) +
  theme(
        legend.title=element_blank(),
        legend.text=element_text(size=24),
        legend.key.size = unit(0.02, 'lines'),
        legend.justification = c(0, 0), 
        legend.position = c(0, 0),
        legend.background = element_rect(fill=NA, colour='black',size=.6, linetype="solid")) +
  geom_text_repel(data=all_results[-log10(all_results$`P-value`)>3 &  all_results$log2FC>1.6,], aes(label=`GeneSymbol`), xlim= c(2, NA) ) + 
  geom_text_repel(data=all_results[-log10(all_results$`P-value`)>3 &  all_results$log2FC< -1.6,], aes(label=`GeneSymbol`), xlim=c(NA,-2) ) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
  
ggsave(all_volc, file = "plots/Volcano Plots for CONT vs. LPS+PYR LA and FC.pdf",
       height=9,width=14,units="in")		  
ggsave(all_volc, file = "plots/Volcano Plots for CONT vs. LPS+PYR LA and FC.tiff",
       height=9,width=14,units="in",dpi=600)		  
##### Venn diagram of overlap between regions with injury model
library(VennDiagram)
FDR_Col = grep("adj.P",colnames(mergedStats),value=T)
SigIDs = sapply(FDR_Col, function(x) mergedStats[mergedStats[,x]<0.10,'transcriptCluster'])

## sig FC
LA_FC_Injury_Venn = venn.diagram(x = SigIDs,
                           category.names = c('FC','LA'),
                           filename = NULL,
                           fill = c(3,4),
                           alpha=0.4,
                           main = NULL,
                           main.cex = 5,
                           cat.just=list(c(0.9,1.5) , c(0.5,0.5) ), 
                           cex=5, cat.cex=4)
tiff('plots/Comparison_Overlap_BetweenRegions_Venn_ID_FDR10.tiff',res=600,height=8,width=8,units='in')
grid.draw(LA_FC_Injury_Venn)
dev.off()
pdf('plots/Comparison_Overlap_BetweenRegions_Venn_ID_FDR10.pdf',height=8,width=8)
grid.draw(LA_FC_Injury_Venn)
dev.off()


#####  Check correlation between t-statistic and signal value
pdf('plots/exploratory_meanSignal_vs_DE_t-value.pdf',height=11,width=8.5)
plot(x=abs(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-t`), y=merged_SstRma_disease$`meanSignal_FC-Control`, xlab="FC: Absolute T-value", ylab='mean signal (Control group)', main = 'FC Injury Effect vs. Control mean expression')
plot( abs(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_FC-Injury`, xlab="FC: Absolute T-value", ylab='mean signal (Injury group)', main = 'FC Injury Effect vs. Injury mean expression' )
plot( abs(merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_LA-Control`, xlab="LA: Absolute T-value", ylab='mean signal (Control group)', main = 'LA Injury Effect vs. Control mean expression')
plot( abs(merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_LA-Injury`, xlab="LA: Absolute T-value", ylab='mean signal (Injury group)', main = 'LA Injury Effect vs. Injury mean expression')
dev.off()

cor(abs(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_FC-Control`)
cor(abs(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_FC-Injury`)
cor(abs(merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_LA-Control`)
cor(abs(merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-t`), merged_SstRma_disease$`meanSignal_LA-Injury`)
