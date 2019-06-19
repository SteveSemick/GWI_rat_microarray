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
lengths(gList)
###
#FC_compareKegg = compareCluster(gList, fun = "enrichKEGG",organism = "rno", 
#                             universe = univ, qvalueCutoff = 0.2, pvalueCutoff = 0.05)
FC_compareGoMf = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
FC_compareGoBp = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Rn.eg.db,qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#FC_compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Rn.eg.db, #qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#compareDO   = compareCluster(gList, fun = "enrichDO", universe = univ,	qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)

FC_compareGo = lapply(list(compareGoMf=FC_compareGoMf, compareGoBp=FC_compareGoBp),simplify)

## Export results
save(FC_compareGo, FC_compareGoMf,FC_compareGoBp, file='rda/FC_GO_Results.rda')
openxlsx::write.xlsx(lapply(FC_compareGo, summary), file = 'csvs/GO_FC_Injury_SST-RMA_R_1e-2_split.xlsx')

## up and down stats - LA
sigStats = merged_SstRma_disease[which(merged_SstRma_disease[,'LA_CTRL_v_PYR_LPS-P.Value'] < 1e-2), c("LA_CTRL_v_PYR_LPS-logFC", "LA_CTRL_v_PYR_LPS-adj.P.Val", "EntrezID")]
sigStats$Sign = sign(sigStats$`LA_CTRL_v_PYR_LPS-logFC`)
sigStats = sigStats[!is.na(sigStats$EntrezID),]
sigStats = sigStats[!duplicated(sigStats[,c("EntrezID")]),]
sigStats =sigStats[-grep(";",sigStats$EntrezID),]
sigStats$EntrezID = as.character(sigStats$EntrezID)
## by sign
gList = split(sigStats$EntrezID, sigStats$Sign)
lengths(gList)
###
LA_compareKegg = compareCluster(gList, fun = "enrichKEGG",organism = "rno", 
                                universe = univ, qvalueCutoff = 0.2, pvalueCutoff = 0.05)
LA_compareGoMf = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
LA_compareGoBp = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Rn.eg.db,qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
LA_compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)

LA_compareGo = lapply(list(compareGoMf=LA_compareGoMf, compareGoBp=LA_compareGoBp, compareGoCc=LA_compareGoCc),simplify)

## Export results
save(LA_compareGo,LA_compareKegg, LA_compareGoMf,LA_compareGoBp, LA_compareGoCc, file='rda/LA_GO_Results.rda')
openxlsx::write.xlsx(lapply(c(KEGG=LA_compareKegg, LA_compareGo), summary), file = 'csvs/GO_LA_Injury_SST-RMA_R_1e-2_split.xlsx')
