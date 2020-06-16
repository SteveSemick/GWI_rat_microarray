# Normalization done with RMA using all 34 samples (both FC and LA; all 4 conditions together).
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/twoBatch_rmaNorm_Data.rda', verbose=T)
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
modFull=model.matrix(~0+pd[,'Group_X_Region']+pd[,'Batch'])
colnames(modFull)= c(levels(pd$Group_X_Region),"Batch2")

#b1 = which(pd$Batch=="Batch1")
#modFull=model.matrix(~pd$Batch+pd$brain_region+pd$Group)
#colnames(modFull)[2:6]= c('Batch2','brain_regionLA','GroupCTRL','GroupPYR_LPS','GroupPYR_LPS_RGZ')

tcMap = tcMap[rownames(normExprs),]

fitFull = lmFit(normExprs, modFull)
conMatInjury = makeContrasts(FC.Injury = PYR_LPS_FC - CTRL_FC,
                             #           PYR_LPS_RGZ_FC-CTRL_FC,
                             #PYR_LPS_RGZ_FC-PYR_LPS_FC,
                             LA.Injury = PYR_LPS_LA-CTRL_LA,
                             # PYR_LPS_RGZ_LA-CTRL_LA,
                             #PYR_LPS_RGZ_LA-PYR_LPS_LA, 
                             levels=fitFull)


fitCon_Injury = contrasts.fit(fitFull,conMatInjury)
fitConEbInjury = eBayes(fitCon_Injury)

results.specific <- decideTests(fitConEbInjury[,c("FC.Injury", "LA.Injury")], method = "global")
summary(results.specific)
Global.Adjusted.P <- fitConEbInjury$p.value
Global.Adjusted.P[] <- p.adjust(Global.Adjusted.P, method="BH")
Global.Adjusted.P = data.frame(Global.Adjusted.P)

## Test models
#Frontal cortex injury models
Full_stats_FC_CTRLvPYR_LPS <- topTable(fitConEbInjury, num = Inf, genelist=tcMap, adjust.method='fdr')
colnames(Full_stats_FC_CTRLvPYR_LPS)[14:17] <- paste0("FC_CTRL_v_PYR_LPS-",colnames(Full_stats_FC_CTRLvPYR_LPS)[14:17])
#Full_stats_FC_CTRLvPYR_LPS$`FC_CTRL_v_PYR_LPS-adj.P.Val` <- Global.Adjusted.P[rownames(Full_stats_FC_CTRLvPYR_LPS), 'FC.Injury']  #Add global-corrected FDR

#Lateral amygdala injury models
Full_stats_LA_CTRLvPYR_LPS <- topTable(fitConEbInjury, num = Inf, coef="LA.Injury", adjust.method='fdr')
colnames(Full_stats_LA_CTRLvPYR_LPS) <- paste0("LA_CTRL_v_PYR_LPS-",colnames(Full_stats_LA_CTRLvPYR_LPS))
#Full_stats_LA_CTRLvPYR_LPS$`LA_CTRL_v_PYR_LPS-adj.P.Val` <- Global.Adjusted.P[rownames(Full_stats_LA_CTRLvPYR_LPS), 'LA.Injury']  #Add global-corrected FDR



#############====== Reorganization and saving =======################
baseRownames = rownames(Full_stats_FC_CTRLvPYR_LPS)
mergedStats=
  cbind(Full_stats_FC_CTRLvPYR_LPS[baseRownames, ],
        Full_stats_LA_CTRLvPYR_LPS[baseRownames, ])
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
save(mergedStats, file='rda/twoBatch_FullModel_mergedStats_RMA_R.rda')
write.csv(mergedStats, file='csvs/twoBatch_FullModel_mergedStats_RMA_R.csv')

