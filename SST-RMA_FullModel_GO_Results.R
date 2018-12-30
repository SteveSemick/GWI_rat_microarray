library(clusterProfiler)
library(DOSE)
library(Reactome)
library(org.Rn.eg.db)
#

load('rda/merged_sstRma_results.rda')
#### 
table(duplicated(merged_SstRma$`Gene Symbol`))
table(grepl(";", merged_SstRma$`Gene Symbol`))

#### GO Analaysis
merged_SstRma_reduced = merged_SstRma[!duplicated(merged_SstRma$`Gene Symbol`) | duplicated(merged_SstRma$`Gene Symbol`,fromLast=T),]
merged_SstRma_reduced = merged_SstRma_reduced[!grepl(";", merged_SstRma_reduced$`Gene Symbol`), ]
library(org.Rn.eg.db)
rn <- org.Rn.eg.db

ann =select(rn, 
            keys = merged_SstRma_reduced$`Gene Symbol`,
            columns = c("ENTREZID"),
            keytype = "SYMBOL")
###
ann = ann[!(duplicated(ann$SYMBOL)|duplicated(ann$SYMBOL,fromLast=T)),]
merged_SstRma_reduced = merged_SstRma_reduced[merged_SstRma_reduced$`Gene Symbol`%in%ann$SYMBOL,]
merged_SstRma_reduced$ENTREZ = ann[match(merged_SstRma_reduced$`Gene Symbol`, ann$SYMBOL),'ENTREZID']
merged_SstRma_reduced = merged_SstRma_reduced[!is.na(merged_SstRma_reduced$`ENTREZ`), ]
###
univ = merged_SstRma_reduced$ENTREZ
univ = as.character(unique(univ[!is.na(univ)]))

Pval_Col = grep("\\.P-val",colnames(merged_SstRma_reduced),value=T)
statList = sapply(Pval_Col, function(x) merged_SstRma_reduced[merged_SstRma_reduced[,x]<5e-3,'ENTREZ'])

library(clusterProfiler)
compareKegg = compareCluster(statList, fun = "enrichKEGG",organism = "rno", universe = univ, qvalueCutoff = 0.2, pvalueCutoff = 0.05)
compareGoMf = compareCluster(statList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(statList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Rn.eg.db,qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(statList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareDO   = compareCluster(statList, fun = "enrichDO", universe = univ,	qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareReact= compareCluster(statList, fun="enrichPathway", universe = univ, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
