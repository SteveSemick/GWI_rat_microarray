## 03 - Testing for differential expression between groups
# Normalization done with SST-RMA using all 16 samples (both FC and LA; all three conditions together).
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/normalized_data_n16_sst-rma.rda')
library(jaffelab)
library(ggplot2)
library(limma)

## Annotation
annotation_Clariom_S_Rat_2018<-read.csv('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/Clariom_S_Rat.r1.na36.rn6.a1.transcript.csv/Clariom_S_Rat.r1.na36.rn6.a1.transcript.csv',colClasses='character',comment.char='#')
a  = strsplit(annotation_Clariom_S_Rat_2018$gene_assignment, split=" /// ")
accession = sapply(a, function(x) paste( jaffelab::ss(x," // ", 1),collapse='; ') )
geneSymbol = sapply(a, function(x) paste( unique(jaffelab::ss(x," // ", 2)) ,collapse='; ') )
gene_description = sapply(a, function(x) paste( unique(jaffelab::ss(x," // ", 3)),collapse='; ') )
EntrezID = gsub("^; ", "", gsub("; $","", sapply(a, function(x) paste( gsub("---","",unique(jaffelab::ss(x," // ", 5))) ,collapse='; ' ) ) ) )
EntrezID[EntrezID==""] <- NA
tcMap = data.frame()
#Clean up strings
EntrezID = EntrezID)



library(clariomsrattranscriptcluster.db)
keytypes(clariomsrattranscriptcluster.db)
select(clariomsrattranscriptcluster.db)
columns(clariomsrattranscriptcluster.db)
probeMap = select(clariomsrattranscriptcluster.db, keys=rownames(sst_rma_signal), columns = c('ENTREZID','GENENAME','SYMBOL' ))

select(clariomsrattranscriptcluster.db, keys='TC0100000006.rn.2', columns = keytypes(clariomsrattranscriptcluster.db))
############
collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
}
# Example:
collapser(ae.annots[ae.annots$PROBEID == dup.ids[1], "SYMBOL"])

# Redefinition of ae.annots
ae.annots <- AnnotationDbi::select(
  x       = moex10sttranscriptcluster.db,
  keys    = rownames(AEset.norm),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  group_by(PROBEID) %>%
  summarise_each(funs(collapser)) %>%
  ungroup

#############
table(rownames(sst_rma_signal)%in%keys(clariomsrattranscriptcluster.db) )
#ID	Gene Symbol	Description	Chromosome	Strand	Start	Stop

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

### Compare to sstRMA
load('TAC/SST-RMA_BothRegions/MergedStats.rda',verbose=T)
