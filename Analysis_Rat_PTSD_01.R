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
## FC 
modFC=model.matrix(~0+pd[pd$brain_region=="FC",'Group'])
colnames(modFC)= c('GroupA', 'GroupB', 'GroupC')
exprsFC=exprs[,which(pd$brain_region=="FC")]

fitFC = lmFit(exprsFC, modFC)
conMat_FC = makeContrasts(GroupA-GroupB,GroupA-GroupC,GroupB-GroupC,levels=fitFC)
fitCon_FC = contrasts.fit(fitFC,conMat_FC)
fitConEb_FC = eBayes(fitCon_FC)

topTable(fitConEb_FC, genelist=probeMap)
topTable(fitConEb_FC, coef="GroupA - GroupB", genelist=probeMap)
topTable(fitConEb_FC, coef="GroupA - GroupC", genelist=probeMap)
topTable(fitConEb_FC, coef="GroupB - GroupC", genelist=probeMap)

### Ordinal model
modFC_ord =model.matrix(~ as.numeric(factor(pd[pd$brain_region=="FC",'Group']) ) )
fitFC_ord = lmFit(exprsFC, modFC_ord)
fitEb_FC_ord = eBayes(fitFC_ord)
res = topTable(fitEb_FC_ord,genelist=probeMap, n=Inf)
table(res$`adj.P.Val`<0.05)

## LA 
modLA=model.matrix(~0+pd[pd$brain_region=="LA",'Group'])
colnames(modLA)= c('GroupA', 'GroupB', 'GroupC')
exprsLA=exprs[,which(pd$brain_region=="LA")]

fitLA = lmFit(exprsLA, modLA)
conMat_LA = makeContrasts(GroupA-GroupB,GroupA-GroupC,GroupB-GroupC,levels=fitLA)
fitCon_LA = contrasts.fit(fitLA,conMat_LA)
fitConEb_LA = eBayes(fitCon_LA)

topTable(fitConEb_LA, genelist=probeMap)[1,]
topTable(fitConEb_LA, coef="GroupA - GroupB",genelist=probeMap)
topTable(fitConEb_LA, coef="GroupA - GroupC",genelist=probeMap)
topTable(fitConEb_LA, coef="GroupB - GroupC",genelist=probeMap)

### Ordinal model
modLA_ord =model.matrix(~ as.numeric(factor(pd[pd$brain_region=="LA",'Group']) ) )
fitLA_ord = lmFit(exprsLA, modLA_ord)
fitEb_LA_ord = eBayes(fitLA_ord)
res = topTable(fitEb_LA_ord,genelist=probeMap, n=Inf)
table(res$`adj.P.Val`<0.05)

### Integrating both brain regions
mod=model.matrix(~0+pd$Group+pd$brain_region)
colnames(mod)= c('GroupA', 'GroupB', 'GroupC', "LA")

fit = lmFit(exprs, mod)
conMat = makeContrasts(GroupA-GroupB,GroupA-GroupC,GroupB-GroupC,levels=fit)
fitCon = contrasts.fit(fit,conMat)
fitConEb = eBayes(fitCon)

res = topTable(fitConEb, genelist=probeMap, n=Inf)
topTable(fitConEb_LA, coef="GroupA - GroupB",genelist=probeMap)
topTable(fitConEb_LA, coef="GroupA - GroupC",genelist=probeMap)
topTable(fitConEb_FC, coef="GroupB - GroupC",genelist=probeMap)

###### Gene set annotation


###
### Create boxplots
pdf('plots/SupplementalFigure_best_ERC_subset_DMC_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (cpg_i in ERC_sigCpG) {
  
  #Change column name for ggplot2 to work
  ii=match(cpg_i, regionSpecific_mergedStats$Name)
  
  pc_genes=unique(unlist(strsplit(regionSpecific_mergedStats[ii,'within10kb_geneSymbol_gencode_hg38'], ";")))
  pc_genes=pc_genes[pc_genes%in%protein_coding_genes_gcV25]
  
  fixName= paste(pc_genes, collapse = ", " )
  custom_title = paste0( 
    "ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
    "\n ",fixName ) #custom title
  
  a = ggplot(dat, aes_string(x = 'Dx', y = cpg_i, fill='Dx')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_jitter(aes(col=`Dx`),width=0.3,size=2 ) + 
    labs(y=paste0(cpg_i, "\nDNAm level"),x="Diagnosis", title = custom_title) + 
    #		scale_colour_brewer(palette = "Set1") +
    #		scale_fill_brewer(palette = "Set1") + 
    scale_fill_manual(values=c("black","#ab1323" ) ) + 
    scale_colour_manual(values=c("black","#ab1323" ) ) + 
    theme(legend.position='none') +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black')) 
  
  print(a)			
}
dev.off()