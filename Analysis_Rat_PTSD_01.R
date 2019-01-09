# Script to import and analyze PTSD Rat Clarion S arrays
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
cel_path = "C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/RatClariomS_16"


##
library(limma)
library(GEOquery)
library(clariomsrattranscriptcluster.db)
## Import .cel files
library(oligo)
library(ggplot2)
library(jaffelab)
library(plyr)
library(tidyr)
theme_set(theme_bw(base_size=40) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))
celFiles = list.files(cel_path,full.names=TRUE,pattern =".CEL")
rawDat = read.celfiles(celFiles)

## Create phenotype table
pd = phenoData(rawDat)
pd = as(pd, "data.frame")
pd$sampleNames = rownames(pd)
pd$brain_region = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",1)

pd$sampleNum = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",2)
pd$sampleNum = as.numeric(gsub("\\..*","",pd$sampleNum))

pd$Group = ifelse(pd$sampleNum %in% c(1,2,3), "GroupA", ifelse(pd$sampleNum %in% c(4,5,6), "GroupB", "GroupC") ) # assigning

## QC on raw data
#perfectMatch = pm(rawDat)
#log_perfectMatch=log2(perfectMatch)+1
#rownames(log_perfectMatch) = 1:nrow(log_perfectMatch)
#dat = cbind(pd,t(log_perfectMatch))
#dat2=gather(dat,key="log_perfectMatch_intensity","value",-index,-sampleNames,-brain_region,-Group,-sampleNum)

## RMA Normalization
library(affycoretools)

data.rma = rma(rawDat)
data.rma <- annotateEset(data.rma, pd.clariom.s.rat)
data.rma = getMainProbes(data.rma) #drop probes not needed

## add annotation using clariom.s.rat package
probeMap = fData(data.rma)

## Extract probe intensities normalized
exprs = exprs(data.rma)

## PCA
pca = prcomp(t(exprs),scale.=TRUE)
varPCs=getPcaVars(pca)
exprsPCs=pca$x
pd=cbind(pd,exprsPCs[,1:10])
pd$PC_VarExplained = varPCs

## PCA plots
PC1_2_region = ggplot(data=pd, aes(x=PC1,y=PC2,col=brain_region ) ) + 
               geom_point(size=5) +	 
               scale_colour_brewer(palette = "Dark2") +
               scale_fill_brewer(palette = "Dark2") + 
               labs(x=paste0("PC1"," (", as.character(signif(pd$PC_VarExplained[1],3)),"%)" ),
                    y=paste0("PC2"," (", as.character(signif(pd$PC_VarExplained[2],3)),"%)" )) +
  theme(legend.position = c(.9,.15), 
        legend.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2")

PC1_2_group = ggplot(data=pd, aes(x=PC1,y=PC2,col=Group ) ) + 
  geom_point(size=5) +	 
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") + 
  labs(x=paste0("PC1"," (", as.character(signif(pd$PC_VarExplained[1],3)),"%)" ),
       y=paste0("PC2"," (", as.character(signif(pd$PC_VarExplained[2],3)),"%)" )) +
  theme(legend.position = c(.85,.15), 
        legend.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Set1")

pdf('plots/PCA_RMA_Normalization.pdf', height=12,width=12)
PC1_2_region
PC1_2_group
dev.off()

### Some stats for PCA
summary(lm(pd$PC1~ as.numeric(as.factor(pd$Group))+pd$brain_region ))

## annotation
#ClariomSAnn = data.table::fread('Clariom_S_Rat.r1.na36.rn6.a1.transcript.csv/Clariom_S_Rat.r1.na36.rn6.a1.transcript.csv')
#probeMap = ClariomSAnn[match(ClariomSAnn$probeset_id, rownames(exprs)),c('probeset_id','seqname','start','stop','gene_assignment')]

#info=ls("package:clariomsrattranscriptcluster.db")
#library(clariomsrattranscriptcluster.db)

#probes=row.names(exprs)
#Symbols = unlist(mget(probes, clariomsrattranscriptclusterSYMBOL, ifnotfound=NA))
#Entrez_IDs = unlist(mget(probes, clariomsrattranscriptclusterENTREZID, ifnotfound=NA))

#probeMap=select(clariomsrattranscriptcluster.db, keys=row.names(exprs), columns=c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","PFAM"), keytype="PROBEID")
#duplicateProbes = probeMap[duplicated(probeMap$PROBEID),'PROBEID']

## Drop probes that map to multiple genes
#probeMap = probeMap[!probeMap$PROBEID%in%duplicateProbes,]
#exprs=exprs[rownames(exprs)%in%probeMap$PROBEID,]
## Modeling with limma
library(limma)

######## Fully partitioned dataset analysis ############3
## FC 
modFC=model.matrix(~0+pd[pd$brain_region=="FC",'Group'])
colnames(modFC)= c('GroupA', 'GroupB', 'GroupC')
exprsFC=exprs[,which(pd$brain_region=="FC")]

fitFC = lmFit(exprsFC, modFC)
conMat_FC = makeContrasts(GroupA-GroupB,GroupA-GroupC,GroupB-GroupC,levels=fitFC)
fitCon_FC = contrasts.fit(fitFC,conMat_FC)
fitConEb_FC = eBayes(fitCon_FC)

## Test models
stats_FC_ALL <- topTable(fitConEb_FC, num = Inf, genelist=probeMap)[,-c(5,6,7)]
colnames(stats_FC_ALL) <- paste0("FC_ALL_",colnames(stats_FC_ALL))

stats_FC_AvB <- topTable(fitConEb_FC, num = Inf, coef="GroupA - GroupB")
colnames(stats_FC_AvB) <- paste0("FC_AvB_",colnames(stats_FC_AvB))

stats_FC_AvC <- topTable(fitConEb_FC, num = Inf, coef="GroupA - GroupC")
colnames(stats_FC_AvC) <- paste0("FC_AvC_",colnames(stats_FC_AvC))

stats_FC_BvC <- topTable(fitConEb_FC, num = Inf, coef="GroupB - GroupC")
colnames(stats_FC_BvC) <- paste0("FC_BvC_",colnames(stats_FC_BvC))

## LA 
modLA=model.matrix(~0+pd[pd$brain_region=="LA",'Group'])
colnames(modLA)= c('GroupA', 'GroupB', 'GroupC')
exprsLA=exprs[,which(pd$brain_region=="LA")]

fitLA = lmFit(exprsLA, modLA)
conMat_LA = makeContrasts(GroupA-GroupB,GroupA-GroupC,GroupB-GroupC,levels=fitLA)
fitCon_LA = contrasts.fit(fitLA,conMat_LA)
fitConEb_LA = eBayes(fitCon_LA)

## Test models
stats_LA_ALL <- topTable(fitConEb_LA, num = Inf)[,-(1:3)]
colnames(stats_LA_ALL) <- paste0("LA_ALL_",colnames(stats_LA_ALL))

stats_LA_AvB <- topTable(fitConEb_LA, num = Inf, coef="GroupA - GroupB")
colnames(stats_LA_AvB) <- paste0("LA_AvB_",colnames(stats_LA_AvB))

stats_LA_AvC <- topTable(fitConEb_LA, num = Inf, coef="GroupA - GroupC")
colnames(stats_LA_AvC) <- paste0("LA_AvC_",colnames(stats_LA_AvC))

stats_LA_BvC <- topTable(fitConEb_LA, num = Inf, coef="GroupB - GroupC")
colnames(stats_LA_BvC) <- paste0("LA_BvC_",colnames(stats_LA_BvC))


############## Merge all these tests ############
baseRownames = rownames(stats_FC_ALL)
mergedStats=
  cbind(stats_FC_ALL[baseRownames, ],
        stats_FC_AvB[baseRownames, ],
        stats_FC_AvC[baseRownames, ],
        stats_FC_BvC[baseRownames, ],
        stats_LA_ALL[baseRownames, ],
        stats_LA_AvB[baseRownames, ],
        stats_LA_AvC[baseRownames, ],
        stats_LA_BvC[baseRownames, ])
sapply(mergedStats[,grep("adj.P.Val",colnames(mergedStats)) ], function(x) sum(x<0.25,na.rm=T) )

### Reorder columns
col_order = c(colnames(mergedStats)[1:4], 
              grep("_adj.P.Val",colnames(mergedStats),value=T),	  
              grep("_logFC",colnames(mergedStats),value=T),
              grep("_P.Value",colnames(mergedStats),value=T),
              grep("_t$",colnames(mergedStats),value=T,fixed=F),
              grep("_AveExpr$",colnames(mergedStats),value=T,fixed=F),
              grep("_B$",colnames(mergedStats),value=T,fixed=F),
              grep("_F$",colnames(mergedStats),value=T,fixed=F))
mergedStats=mergedStats[,col_order]			
save(mergedStats, file='rda/mergedStats_RMA_R.rda')
write.csv(mergedStats, file='csvs/mergedStats_RMA_R.csv')


######## Contrasted dataset analysis ############3
pd$Group_X_Region = paste0(pd$Group, "_", pd$brain_region)
pd$Group_X_Region = factor(pd$Group_X_Region)
modFull=model.matrix(~0+pd[,'Group_X_Region'])
colnames(modFull)= levels(pd$Group_X_Region)

fitFull = lmFit(exprs, modFull)
conMatFull = makeContrasts(GroupA_FC-GroupB_FC,GroupA_FC-GroupC_FC,GroupB_FC-GroupC_FC,GroupA_LA-GroupB_LA,GroupA_LA-GroupC_LA,GroupB_LA-GroupC_LA, levels=fitFull)
fitCon_Full = contrasts.fit(fitFull,conMatFull)
fitConEb = eBayes(fitCon_Full)

## Test models
Full_stats_FC_AvB <- topTable(fitConEb, num = Inf, coef="GroupA_FC - GroupB_FC", genelist=probeMap)
colnames(Full_stats_FC_AvB) <- paste0("FC_AvB_",colnames(Full_stats_FC_AvB))

Full_stats_FC_AvC <- topTable(fitConEb, num = Inf, coef="GroupA_FC - GroupC_FC")
colnames(Full_stats_FC_AvC) <- paste0("FC_AvC_",colnames(Full_stats_FC_AvC))

Full_stats_FC_BvC <- topTable(fitConEb, num = Inf, coef="GroupB_FC - GroupC_FC")
colnames(Full_stats_FC_BvC) <- paste0("FC_BvC_",colnames(Full_stats_FC_BvC))

Full_stats_LA_AvB <- topTable(fitConEb, num = Inf, coef="GroupA_LA - GroupB_LA")
colnames(Full_stats_LA_AvB) <- paste0("LA_AvB_",colnames(Full_stats_LA_AvB))

Full_stats_LA_AvC <- topTable(fitConEb, num = Inf, coef="GroupA_LA - GroupC_LA")
colnames(Full_stats_LA_AvC) <- paste0("LA_AvC_",colnames(Full_stats_LA_AvC))

Full_stats_LA_BvC <- topTable(fitConEb, num = Inf, coef="GroupB_LA - GroupC_LA")
colnames(Full_stats_LA_BvC) <- paste0("LA_BvC_",colnames(Full_stats_LA_BvC))

############## Merge all these tests ############
baseRownames = rownames(Full_stats_FC_AvB)
mergedStats=
  cbind(Full_stats_FC_AvB[baseRownames, ],
        Full_stats_FC_AvC[baseRownames, ],
        Full_stats_FC_BvC[baseRownames, ],
        Full_stats_LA_AvB[baseRownames, ],
        Full_stats_LA_AvC[baseRownames, ],
        Full_stats_LA_BvC[baseRownames, ])
sapply(mergedStats[,grep("adj.P.Val",colnames(mergedStats)) ], function(x) sum(x<0.10,na.rm=T) )

### Reorder columns
col_order = c(colnames(mergedStats)[1:4], 
              grep("_adj.P.Val",colnames(mergedStats),value=T),	  
              grep("_logFC",colnames(mergedStats),value=T),
              grep("_P.Value",colnames(mergedStats),value=T),
              grep("_t$",colnames(mergedStats),value=T,fixed=F),
              grep("_AveExpr$",colnames(mergedStats),value=T,fixed=F),
              grep("_B$",colnames(mergedStats),value=T,fixed=F),
              grep("_F$",colnames(mergedStats),value=T,fixed=F))
mergedStats=mergedStats[,col_order]			
save(mergedStats, file='rda/FullModel_mergedStats_RMA_R.rda')
write.csv(mergedStats, file='csvs/FullModel_mergedStats_RMA_R.csv')

###### Gene set annotation

