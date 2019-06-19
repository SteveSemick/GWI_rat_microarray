## 04b - Boxplots for differentially expressed genes between 
# Plotting done with the results of differential expression in R between SST-RMA normalized signals using all 16 samples (both FC and LA; all three conditions together).

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

##

dat = cbind(pd, t(sst_rma_signal)  )
dat = dat[dat$Group%in% c("CTRL","PYR_LPS"),]

merged_SstRma_disease = merged_SstRma_disease[order(merged_SstRma_disease$Region_MinP),]
injury_genes_plot = as.character ( merged_SstRma_disease[merged_SstRma_disease$LA_sig|merged_SstRma_disease$FC_sig,'transcriptCluster']) 

pdf('plots/InjuryGenes_boxplots.pdf',height=10,width=14,useDingbats=FALSE)
for (probe_i in injury_genes_plot) {
  
  #Change column name for ggplot2 to work
  ii=match(probe_i, merged_SstRma_disease$transcriptCluster)
  
  #  custom_title = paste0( 
  #    "ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
  #    "\n ",fixName ) #custom title
  custom_title=paste0(merged_SstRma_disease[ii,'Symbol'], "\nMinP = ", as.character(signif(merged_SstRma_disease[ii,'Region_MinP'],3)) )
  
  a = ggplot(dat, aes_string(x = 'brain_region', y = probe_i, fill='Group')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_point(aes(col=`Group`),size=5,position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
    labs(y=paste0(probe_i, "\nNormalized Signal"),x='Brain Region', title = custom_title) + 
    theme(legend.position='bottom') +
  #  scale_colour_brewer(palette = "Dark2") +
 #   scale_fill_brewer(palette = "Dark2") + 
    scale_fill_manual(values=c("black","#ab1323" ) ) + 
    scale_colour_manual(values=c("black","#ab1323" ) ) + 
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_text(size=40,colour='black'),
          legend.position='none')
  print(a)			
}
dev.off()

#Create specific tiff plot

tiff('plots/Fig1b_Chrm5_boxplots.tiff',height=10,width=14, units='in', res=600)
probe_i='TC0300003893.rn.2'
  #Change column name for ggplot2 to work
  ii=match(probe_i, merged_SstRma_disease$transcriptCluster)
  
  custom_title=paste0(merged_SstRma_disease[ii,'Symbol'] )
  
  a = ggplot(dat, aes_string(x = 'brain_region', y = probe_i, fill='Group')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_point(aes(col=`Group`),size=5,position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
    labs(y=paste0(probe_i, "\nNormalized Signal"),x='Brain Region', title = custom_title) + 
    theme(legend.position='bottom') +
    #  scale_colour_brewer(palette = "Dark2") +
    #   scale_fill_brewer(palette = "Dark2") + 
    scale_fill_manual(values=c("black","#ab1323" ) ) + 
    scale_colour_manual(values=c("black","#ab1323" ) ) + 
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_text(size=40,colour='black'),
          legend.position='none')
  print(a)			
dev.off()