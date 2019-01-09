RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/SST-RMA_BothRegions',pattern = '.txt',full.names = T)

res = lapply(RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(RMA_Res))

### Extract Normalized Signal
rma_signal = res[[1]]
rownames(rma_signal) <- rma_signal$ID
rma_signal = rma_signal[, grep("Signal", colnames(rma_signal)) ]

pd=data.frame(sampleNames=colnames(rma_signal),stringsAsFactors=FALSE)
pd$brain_region = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",1)

pd$sampleNum = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",2)
pd$sampleNum = as.numeric(gsub("\\..*","",pd$sampleNum))

pd$Group = ifelse(pd$sampleNum %in% c(1,2,3), "GroupA", ifelse(pd$sampleNum %in% c(4,5,6), "GroupB", "GroupC") ) # assigning


injury_genes <- read.csv('csvs/Injury_genes_of_interest.csv')
injury_genes_plot = injury_genes$ID

pd$Group = plyr::revalue(pd$Group, c("GroupA"="Control", "GroupB"="Injury", "GroupC"="Injury + Trt"))
dat = cbind(pd, t(rma_signal)  )

pdf('plots/InjuryGenes_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (probe_i in injury_genes_plot) {
  
  #Change column name for ggplot2 to work
  ii=match(probe_i, injury_genes$ID)
  
#  custom_title = paste0( 
#    "ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
#    "\n ",fixName ) #custom title
 custom_title=injury_genes[ii,'Gene.Symbol']
 
  a = ggplot(dat, aes_string(x = 'Group', y = probe_i, fill='brain_region')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_point(aes(col=`brain_region`),size=5,position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
    labs(y=paste0(probe_i, " Log2(Signal)"),x='Group', title = custom_title) + 
    theme(legend.position='bottom') +
    		scale_colour_brewer(palette = "Dark2") +
    		scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
  print(a)			
}
dev.off()

####
LA_rescued_genes <- read.csv('csvs/LA_rescued_genes_of_interest.csv')
LA_genes_plot = LA_rescued_genes$ID

pdf('plots/LateralAmygdala_rescued_Genes_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (probe_i in LA_genes_plot) {
  
  #Change column name for ggplot2 to work
  ii=match(probe_i, LA_rescued_genes$ID)
  
  #  custom_title = paste0( 
  #    "ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
  #    "\n ",fixName ) #custom title
  custom_title=LA_rescued_genes[ii,'Gene.Symbol']
  
  a = ggplot(dat, aes_string(x = 'Group', y = probe_i, fill='brain_region')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_point(aes(col=`brain_region`),size=5,position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
    labs(y=paste0(probe_i, " Log2(Signal)"),x='Group', title = custom_title) + 
    theme(legend.position='bottom') +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
  print(a)			
}
dev.off()

####
FC_rescued_genes <- read.csv('csvs/FC_rescued_genes_of_interest.csv')
FC_genes_plot = FC_rescued_genes$ID

pdf('plots/FrontalCortex_rescued_Genes_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (probe_i in FC_genes_plot) {
  
  #Change column name for ggplot2 to work
  ii=match(probe_i, FC_rescued_genes$ID)
  
  #  custom_title = paste0( 
  #    "ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
  #    "\n ",fixName ) #custom title
  custom_title=FC_rescued_genes[ii,'Gene.Symbol']
  
  a = ggplot(dat, aes_string(x = 'Group', y = probe_i, fill='brain_region')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_point(aes(col=`brain_region`),size=5,position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
    labs(y=paste0(probe_i, " Log2(Signal)"),x='Group', title = custom_title) + 
    theme(legend.position='bottom') +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
  print(a)			
}
dev.off()