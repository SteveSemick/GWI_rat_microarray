load('rda/merged_sstRma_results.rda',verbose=T)
library(dplyr)
library(ggplot2)
library('ggrepel')
theme_set(theme_bw(base_size=24) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"				 )) 




colnames(merged_SstRma)
merged_SstRma_disease = merged_SstRma[,c(1, 12:19, grep("AvB", colnames(merged_SstRma) ) ) ]

merged_SstRma_disease

merged_SstRma_disease$LA_AvB_log2FC = merged_SstRma_disease[,'LA_AvB.LA: A Avg (log2)'] - merged_SstRma_disease[,'LA_AvB.LA: B Avg (log2)']
merged_SstRma_disease$FC_AvB_log2FC = merged_SstRma_disease[,'FC_AvB.FC: A Avg (log2)'] - merged_SstRma_disease[,'FC_AvB.FC: B Avg (log2)']

merged_SstRma_disease$LA_AvB_Bonf_P = p.adjust(merged_SstRma_disease$`LA_AvB.P-val`,method='bonferroni')
merged_SstRma_disease$FC_AvB_Bonf_P = p.adjust(merged_SstRma_disease$`FC_AvB.P-val`,method='bonferroni')

####### Cross-region comparisons of injury  
merged_SstRma_disease = dplyr::mutate(as.data.frame(merged_SstRma_disease), LA_sig=ifelse(merged_SstRma_disease$`LA_AvB.FDR P-val`< 0.10, TRUE, FALSE))
merged_SstRma_disease = dplyr::mutate(as.data.frame(merged_SstRma_disease), FC_sig=ifelse(merged_SstRma_disease$`FC_AvB.FDR P-val`< 0.10, TRUE, FALSE))
merged_SstRma_disease$Sig = "Not sig"
merged_SstRma_disease$Sig[merged_SstRma_disease$FC_sig] = "FC Sig"
merged_SstRma_disease$Sig[merged_SstRma_disease$LA_sig] = "LA Sig"
merged_SstRma_disease$Sig[merged_SstRma_disease$LA_sig &  merged_SstRma_disease$FC_sig] = "Both Sig"
merged_SstRma_disease$Sig = as.factor(merged_SstRma_disease$Sig)

lfc_comp = ggplot(data=merged_SstRma_disease, aes(x=FC_AvB_log2FC, y= LA_AvB_log2FC, col=Sig ))				 
lfc_comp = lfc_comp + geom_point()+ labs(x="Frontal Cortex log2 fold change", y="Amygdala log2 fold change", title = 'Control vs. Injury Effects') + 
  geom_hline(yintercept=0,colour='red',size=1) + geom_vline(xintercept=0,colour='red',size=1) +
  scale_colour_manual(values=c("red","steelblue", "green3","grey") ) + theme(legend.position='bottom')
genes_of_high_interest = merged_SstRma_disease[(merged_SstRma_disease$FC_AvB_log2FC>1.5 & merged_SstRma_disease$LA_AvB_log2FC>1.5) | (merged_SstRma_disease$FC_AvB_log2FC < -1.5 & merged_SstRma_disease$LA_AvB_log2FC< -1.5),]
lfc_comp = lfc_comp + geom_text_repel(data=genes_of_high_interest, aes(label=`Gene Symbol`))
ggsave(lfc_comp, file = "plots/Control vs Injury effect size scatter 1pt5LFC threshold.pdf",
       height=8.5,width=11,units="in")		   

genes_of_high_interest$Region_MinPval = matrixStats::rowMins(as.matrix(genes_of_high_interest[,c('FC_AvB.P-val', 'LA_AvB.P-val' )]) )
genes_of_high_interest[order(genes_of_high_interest$Region_MinPval),c('Gene Symbol','Region_MinPval')]
write.csv(genes_of_high_interest[order(genes_of_high_interest$Region_MinPval),c('ID','Gene Symbol','Description','Chromosome','Strand','Start','Stop', 'Region_MinPval','FC_AvB.P-val','LA_AvB.P-val','LA_AvB_log2FC', 'FC_AvB_log2FC')],file='csvs/Injury_genes_of_interest.csv',row.names=FALSE)

plot(merged_SstRma_disease$LA_AvB_log2FC, merged_SstRma_disease$FC_AvB_log2FC)
cor.test(merged_SstRma_disease$LA_AvB_log2FC, merged_SstRma_disease$FC_AvB_log2FC)

plot(-log10(merged_SstRma_disease$`LA_AvB.P-val`), -log10(merged_SstRma_disease$`FC_AvB.P-val`) )
cor.test(merged_SstRma_disease$LA_AvB_log2FC, merged_SstRma_disease$FC_AvB_log2FC)

############################### Volcano Facet Plot
LA_max_P_FDR10 = max(merged_SstRma_disease[merged_SstRma_disease$`LA_AvB.FDR P-val`<.10,'LA_AvB.P-val'])
FC_max_P_FDR10 = max(merged_SstRma_disease[merged_SstRma_disease$`FC_AvB.FDR P-val`<.10,'FC_AvB.P-val'])
shared_thresh = min(FC_max_P_FDR10,LA_max_P_FDR10)

LA_results = dplyr::mutate(as.data.frame(merged_SstRma_disease), sig=ifelse(merged_SstRma_disease$`LA_AvB.FDR P-val`< 0.10, "FDR<0.10", "FDR>0.10"))
LA_results <- LA_results[,c("Gene Symbol","LA_AvB.P-val","LA_AvB_log2FC","sig")]
colnames(LA_results) = c("GeneSymbol", "P-value", "log2FC", "Sig")
LA_results$Type = "Lateral Amygdala"

FC_results = dplyr::mutate(as.data.frame(merged_SstRma_disease), sig=ifelse(merged_SstRma_disease$`FC_AvB.FDR P-val`< 0.10, "FDR<0.10", "FDR>0.10"))
FC_results <- FC_results[,c("Gene Symbol","FC_AvB.P-val","FC_AvB_log2FC","sig")]
colnames(FC_results) = c("GeneSymbol", "P-value", "log2FC", "Sig")
FC_results$Type = "Frontal Cortex"

all_results <- rbind(LA_results, FC_results)

## Plotting
all_volc <- ggplot(all_results, aes(log2FC, -log10(`P-value`))) +
  geom_point(aes(col=Sig)) +
  facet_wrap(~ Type, ncol= 2, nrow=1) +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.background = element_rect(fill=NA, size=.5, linetype="solid")) +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"(pvalue)")) + 
  scale_color_manual(values=c("red", "grey")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, round(max(-log10(all_results$`P-value`))+.5) ) ) +
  theme(legend.background = element_rect(colour = "black"),
        legend.title=element_blank()) #no legend title 

#all_volc <- all_volc + geom_text_repel(data=all_results[all_results[,'Sig']=="FDR<0.10",], aes(label=GeneSymbol))

ggsave(all_volc, file = "plots/Volcano Plots for CONT vs. Injury LA and FC.pdf",
       height=8.5,width=11,units="in")		   
