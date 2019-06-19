### 06 
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/normalized_data_n16_sst-rma.rda', verbose=T)
library(jaffelab)
library(ggplot2)
library(limma)

#############====== Annotation =======################

## Parse annotation provided by Affymetrix
annotation_Clariom_S_Rat_2018<-read.csv('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/Clariom_S_Rat.r1.na36.rn6.a1.transcript.csv/Clariom_S_Rat.r1.na36.rn6.a1.transcript.csv',colClasses='character',comment.char='#')
a  = strsplit(annotation_Clariom_S_Rat_2018$gene_assignment, split=" /// ")
accession = sapply(a, function(x) paste( jaffelab::ss(x," // ", 1),collapse='; ') )
geneSymbol = gsub("^; ", "", sapply(a, function(x) paste( gsub("---","",unique(jaffelab::ss(x," // ", 2))) ,collapse='; ' )  ) )
gene_description = sapply(a, function(x) paste( unique(jaffelab::ss(x," // ", 3)),collapse='; ') )
EntrezID = gsub("^; ", "", gsub("; $","", sapply(a, function(x) paste( gsub("---","",unique(jaffelab::ss(x," // ", 5))) ,collapse='; ' ) ) ) )
EntrezID[EntrezID==""] <- NA
#Combine annotation into a useful dataframe
tcMap = data.frame(transcriptCluster=annotation_Clariom_S_Rat_2018$transcript_cluster_id,
                   probeset=annotation_Clariom_S_Rat_2018$probeset_id,
                   Chr=annotation_Clariom_S_Rat_2018$seqname,
                   Strand=annotation_Clariom_S_Rat_2018$strand,
                   Start = annotation_Clariom_S_Rat_2018$start,
                   End = annotation_Clariom_S_Rat_2018$stop,
                   Accession = accession,
                   Symbol=geneSymbol,
                   EntrezID=EntrezID,
                   gene_description=gene_description,
                   locusType = annotation_Clariom_S_Rat_2018$locus.type)
rownames(tcMap) <-tcMap$transcriptCluster

#############====== Differential expression testing =======################
##set up
pd$Group_X_Region = paste0(pd$Group, "_", pd$brain_region)
pd$Group_X_Region = factor(pd$Group_X_Region)
modFull=model.matrix(~0+pd[,'Group_X_Region'])
colnames(modFull)= levels(pd$Group_X_Region)

tcMap = tcMap[rownames(sst_rma_signal),]

fitFull = lmFit(sst_rma_signal, modFull)
# Test only for differences between roglitizone treatment and injury
conMatTreat = makeContrasts(FC.Treat = PYR_LPS_RGZ_FC-PYR_LPS_FC,
                             LA.Treat = PYR_LPS_RGZ_LA-PYR_LPS_LA,
                             levels=fitFull)

fitCon_Treat = contrasts.fit(fitFull,conMatTreat)
fitConEbTreat = eBayes(fitCon_Treat)

### Subset to only genes disrupted with injury
load('rda/injury_DE_res.rda', verbose=T)
LA_tc_of_interest = as.character(merged_SstRma_disease[merged_SstRma_disease$LA_sig, 'transcriptCluster' ])
LA_indicies_of_interest = match(LA_tc_of_interest, rownames(sst_rma_signal) )

FC_tc_of_interest = as.character(merged_SstRma_disease[merged_SstRma_disease$FC_sig, 'transcriptCluster' ])
LA_indicies_of_interest = match(FC_tc_of_interest, rownames(sst_rma_signal) )


results.specific <- decideTests(fitConEbTreat[indicies_of_interest,c("FC.Treat", "LA.Treat")], method = "global")
summary(results.specific)
Global.Adjusted.P <- fitConEbTreat[indicies_of_interest,]$p.value
Global.Adjusted.P[] <- p.adjust(Global.Adjusted.P, method="BH")
Global.Adjusted.P = data.frame(Global.Adjusted.P)

##
#Frontal cortex treat models
FC_Full_stats_Treat <- topTable(fitConEbTreat[indicies_of_interest,], num = Inf, coef="FC.Treat", genelist=tcMap, adjust.method='BH')
colnames(FC_Full_stats_Treat)[12:17] <- paste0("FC_Treat-",colnames(FC_Full_stats_Treat)[12:17])
FC_Full_stats_Treat$`FC_Treat-adj.P.Val` <- Global.Adjusted.P[rownames(FC_Full_stats_Treat), 'FC.Treat']  #Add global-corrected FDR

#Lateral amygdala treat models
LA_Full_stats_Treat <- topTable(fitConEbTreat[indicies_of_interest,], num = Inf, coef="LA.Treat", adjust.method='BH')
colnames(LA_Full_stats_Treat) <- paste0("LA_Treat-",colnames(LA_Full_stats_Treat) )
LA_Full_stats_Treat$`LA_Treat-adj.P.Val` <- Global.Adjusted.P[rownames(LA_Full_stats_Treat), 'LA.Treat']  #Add global-corrected FDR

#############====== Reorganization and saving =======################
baseRownames = rownames(FC_Full_stats_Treat)
mergedStats=
  cbind(FC_Full_stats_Treat[baseRownames, ],
        LA_Full_stats_Treat[baseRownames, ])
sapply(mergedStats[,grep("adj.P.Val",colnames(mergedStats)) ], function(x) sum(x<0.05,na.rm=T) )

### Reorder columns
col_order = c(colnames(mergedStats)[1:11], 
              grep("-adj.P.Val",colnames(mergedStats),value=T),	  
              grep("-logFC",colnames(mergedStats),value=T),
              grep("-P.Value",colnames(mergedStats),value=T),
              grep("-t$",colnames(mergedStats),value=T,fixed=F),
              grep("-AveExpr$",colnames(mergedStats),value=T,fixed=F),
              grep("-B$",colnames(mergedStats),value=T,fixed=F))
mergedStats=mergedStats[,col_order]		
Treat_mergedStats = mergedStats
save(Treat_mergedStats, file='rda/RGZ_Treat_FullModel_mergedStats_sst-RMA_R.rda')
write.csv(Treat_mergedStats, file='csvs/RGZ_Treat_FullModel_mergedStats_sst-RMA_R.csv')
#####
allConditions_Merged = cbind(merged_SstRma_disease[rownames(Treat_mergedStats), ], Treat_mergedStats[,11:23 ])
  

### Scatterplot of effect sizes between injury and treatment
library(jaffelab)
library(ggplot2)
library(ggrepel)
library(limma)
theme_set(theme_bw(base_size=40) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))
fc_lfc_comp = ggplot(data=allConditions_Merged, aes(x=`FC_CTRL_v_PYR_LPS-logFC`, y= `FC_Treat-logFC`))				 
fc_lfc_comp = fc_lfc_comp + 
  geom_point() + 
  labs(x="FC Injury log2 fold change", y="FC Treatment log2 fold change") + 
  geom_hline(yintercept=0,colour='red',size=1) + 
  geom_vline(xintercept=0,colour='red',size=1)+ 
  geom_abline(colour='steelblue',size=1, slope=-1)
# scale_colour_manual(values=c("red","steelblue", "green3","grey") ) + 
#genes_of_high_interest = merged_SstRma_disease[(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC`>1.5 & merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`>1.5) | (merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC` < -1.5 & merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`< -1.5),]
#lfc_comp = lfc_comp + geom_text_repel(data=genes_of_high_interest, aes(label=`Symbol`))

ggsave(fc_lfc_comp, file = "plots/FC Injury vs. Treatment Effect Sizes.pdf",
       height=12,width=12,units="in")		   

ggsave(fc_lfc_comp, file = "plots/FC Injury vs. Treatment Effect Sizes.tiff",
       height=12,width=12,units="in", dpi=600)		  

######### Lateral amygdala fold change effects
la_lfc_comp = ggplot(data=allConditions_Merged, aes(x=`LA_CTRL_v_PYR_LPS-logFC`, y= `LA_Treat-logFC`))				 
la_lfc_comp = la_lfc_comp + 
  geom_point() + 
  labs(x="LA Injury log2 fold change", y="LA Treatment log2 fold change") + 
  geom_hline(yintercept=0,colour='red',size=1) + 
  geom_vline(xintercept=0,colour='red',size=1) + 
  geom_abline(colour='steelblue',size=1, slope=-1)
# scale_colour_manual(values=c("red","steelblue", "green3","grey") ) + 
#genes_of_high_interest = merged_SstRma_disease[(merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC`>1.5 & merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`>1.5) | (merged_SstRma_disease$`FC_CTRL_v_PYR_LPS-logFC` < -1.5 & merged_SstRma_disease$`LA_CTRL_v_PYR_LPS-logFC`< -1.5),]
#lfc_comp = lfc_comp + geom_text_repel(data=genes_of_high_interest, aes(label=`Symbol`))

ggsave(la_lfc_comp, file = "plots/LA Injury vs. Treatment Effect Sizes.pdf",
       height=12,width=12,units="in")		   

ggsave(la_lfc_comp, file = "plots/LA Injury vs. Treatment Effect Sizes.tiff",
       height=12,width=12,units="in", dpi=600)		  
