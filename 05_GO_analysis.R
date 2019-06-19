## 05 - Gene set enrichment analysis
# Gene set enrichment analysis 

library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(org.Rn.eg.db)
#
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/injury_DE_res.rda', verbose=T)
###
univ = merged_SstRma_disease$EntrezID
univ = as.character(unique(univ[!is.na(univ)]))
univ= univ[-grep(";",univ)]#drop TCs w/ multiple entrezIDs
length(univ)
###

## up and down stats - FC
Pval_Col = grep("P.Value",colnames(merged_SstRma_disease),value=T)
sigStats = merged_SstRma_disease[which(merged_SstRma_disease[,'FC_CTRL_v_PYR_LPS-P.Value'] < 1e-2), c("FC_CTRL_v_PYR_LPS-logFC", "FC_CTRL_v_PYR_LPS-adj.P.Val", "EntrezID")]
sigStats$Sign = sign(sigStats$`FC_CTRL_v_PYR_LPS-logFC`)
sigStats = sigStats[!is.na(sigStats$EntrezID),]
sigStats = sigStats[!duplicated(sigStats[,c("EntrezID")]),]
sigStats =sigStats[-grep(";",sigStats$EntrezID),]
sigStats$EntrezID = as.character(sigStats$EntrezID)
## by sign
gList = split(sigStats$EntrezID, sigStats$Sign)
lengths(gLsit)

compareKegg = compareCluster(statList, fun = "enrichKEGG",organism = "rno", universe = univ, qvalueCutoff = 0.2, pvalueCutoff = 0.05)
compareGoMf = compareCluster(statList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(statList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Rn.eg.db,qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(statList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#compareDO   = compareCluster(statList, fun = "enrichDO", universe = univ,	qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#compareReact= compareCluster(statList, fun="enrichPathway", universe = univ, 
#                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)

compareGo = lapply(list(compareGoMf=compareGoMf, compareGoBp=compareGoBp, compareGoCc=compareGoCc),simplify)

## Export results
save(compareKegg,compareGoMf,compareGoBp,compareGoCc, file='rda/GO_Results.rda')
openxlsx::write.xlsx(lapply(c(compareGo,compareKegg=compareKegg), summary),
                     file = 'csvs/GO_Analysis_PTSD_FullModel_SST-RMA_TAC_5e-3_unsplit.xlsx')
