load('rda/merged_sstRma_results.rda',verbose=T)
library(dplyr)
library(ggplot2)
library('ggrepel')
theme_set(theme_bw(base_size=24) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"				 )) 



### Calculating LFC
merged_SstRma$LA_AvB_log2FC = merged_SstRma[,'LA_AvB.LA: A Avg (log2)'] - merged_SstRma[,'LA_AvB.LA: B Avg (log2)']
merged_SstRma$FC_AvB_log2FC = merged_SstRma[,'FC_AvB.FC: A Avg (log2)'] - merged_SstRma[,'FC_AvB.FC: B Avg (log2)']

merged_SstRma$LA_BvC_log2FC = merged_SstRma[,'LA_BvC.LA: B Avg (log2)'] - merged_SstRma[,'LA_BvC.LA: C Avg (log2)']
merged_SstRma$FC_BvC_log2FC = merged_SstRma[,'FC_BvC.FC: B Avg (log2)'] - merged_SstRma[,'FC_BvC.FC: C Avg (log2)']

### Subsetting to genes 
# A-B >0 => B<A => injury reduces genes expression
# B-C <0 => B<C => treatment increases gene expression


merged_SstRma = dplyr::mutate(as.data.frame(merged_SstRma), LA_AvB_sig=ifelse(merged_SstRma$`LA_AvB.P-val`< 1e-3, TRUE, FALSE))
merged_SstRma = dplyr::mutate(as.data.frame(merged_SstRma), FC_AvB_sig=ifelse(merged_SstRma$`FC_AvB.P-val`< 1e-3, TRUE, FALSE))
merged_SstRma = dplyr::mutate(as.data.frame(merged_SstRma), LA_BvC_sig=ifelse(merged_SstRma$`LA_BvC.P-val`< 1e-3, TRUE, FALSE))
merged_SstRma = dplyr::mutate(as.data.frame(merged_SstRma), FC_BvC_sig=ifelse(merged_SstRma$`FC_BvC.P-val`< 1e-3, TRUE, FALSE))

merged_SstRma$FC_Sig = "Not sig"
merged_SstRma$FC_Sig[merged_SstRma$FC_AvB_sig] = "Injury Sig"
merged_SstRma$FC_Sig[merged_SstRma$FC_BvC_sig] = "Treatment Sig"
merged_SstRma$FC_Sig[merged_SstRma$FC_AvB_sig &  merged_SstRma$FC_BvC_sig] = "Injury and Treatment Sig"
merged_SstRma$FC_Sig = factor(merged_SstRma$FC_Sig, c("Injury and Treatment Sig", "Injury Sig", "Treatment Sig", "Not sig") )

lfc_comp_FC = ggplot(data=merged_SstRma, aes(x=FC_AvB_log2FC, y= FC_BvC_log2FC,col=FC_Sig ))				 
lfc_comp_FC = lfc_comp_FC + geom_point()+ labs(x="Control vs. Injury LFC", y="Injury vs. Trt LFC", title = 'Frontal Cortex Effect Sizes') + 
  geom_hline(yintercept=0,colour='red',size=1) + geom_vline(xintercept=0,colour='red',size=1) +
  scale_colour_manual(values=c("red","steelblue", "green3","grey") ) + theme(legend.position='bottom')

genes_of_high_interest_FC = merged_SstRma[(merged_SstRma$FC_AvB_log2FC< -1.5 & merged_SstRma$FC_BvC_log2FC>1) | (merged_SstRma$FC_AvB_log2FC > 1.5 & merged_SstRma$FC_BvC_log2FC< -1),]
lfc_comp_FC = lfc_comp_FC + geom_text_repel(data=genes_of_high_interest_FC, aes(label=`Gene Symbol`))
write.csv(genes_of_high_interest_FC[order(genes_of_high_interest_FC[,c('FC_AvB.P-val')],genes_of_high_interest_FC[,('FC_BvC.P-val')]),c('ID','Gene Symbol','Description','Chromosome','Strand','Start','Stop','FC_AvB.P-val','FC_BvC.P-val','FC_AvB_log2FC', 'FC_BvC_log2FC')],file='csvs/FC_rescued_genes_of_interest.csv',row.names=FALSE)

ggsave(lfc_comp_FC, file = "plots/FC Rescue expression effect size scatter 1pt5LFC threshold.pdf",
       height=8.5,width=11,units="in")		   

############ LATERAL AMYGDALA ##################
merged_SstRma$LA_Sig = "Not sig"
merged_SstRma$LA_Sig[merged_SstRma$LA_AvB_sig] = "Injury Sig"
merged_SstRma$LA_Sig[merged_SstRma$LA_BvC_sig] = "Treatment Sig"
merged_SstRma$LA_Sig[merged_SstRma$LA_AvB_sig &  merged_SstRma$LA_BvC_sig] = "Injury and Treatment Sig"
merged_SstRma$LA_Sig = factor(merged_SstRma$LA_Sig, c("Injury and Treatment Sig", "Injury Sig", "Treatment Sig", "Not sig") )

lfc_comp_LA = ggplot(data=merged_SstRma, aes(x=LA_AvB_log2FC, y= LA_BvC_log2FC,col=LA_Sig ))				 
lfc_comp_LA = lfc_comp_LA + geom_point()+ labs(x="Control vs. Injury LFC", y="Injury vs. Trt LFC", title = 'Lateral Amygdala Effect Sizes') + 
  geom_hline(yintercept=0,colour='red',size=1) + geom_vline(xintercept=0,colour='red',size=1) +
  scale_colour_manual(values=c("red","steelblue", "green3","grey") ) + theme(legend.position='bottom')

genes_of_high_interest_LA = merged_SstRma[(merged_SstRma$LA_AvB_log2FC< -1.5 & merged_SstRma$LA_BvC_log2FC>1) | (merged_SstRma$LA_AvB_log2FC > 1.5 & merged_SstRma$LA_BvC_log2FC< -1),]
lfc_comp_LA = lfc_comp_LA + geom_text_repel(data=genes_of_high_interest_LA, aes(label=`Gene Symbol`))
write.csv(genes_of_high_interest_LA[order(genes_of_high_interest_LA[,c('LA_AvB.P-val')],genes_of_high_interest_LA[,c('LA_BvC.P-val')]),c('ID','Gene Symbol','Description','Chromosome','Strand','Start','Stop','LA_AvB.P-val','LA_BvC.P-val','LA_AvB_log2FC', 'LA_BvC_log2FC')],file='csvs/LA_rescued_genes_of_interest.csv',row.names=FALSE)

ggsave(lfc_comp_LA, file = "plots/LA Rescue expression effect size scatter 1pt5LFC threshold.pdf",
       height=8.5,width=11,units="in")		   