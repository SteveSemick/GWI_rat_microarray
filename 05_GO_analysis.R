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


## Check for overlap between FC and LA Go categories
load('rda/LA_GO_Results.rda')
load('rda/FC_GO_Results.rda',verbose=T)

intersect(summary(LA_compareGoMf)$ID, summary(FC_compareGoMf)$ID)

simple_LA_compareGoBp = simplify(LA_compareGoBp)
simple_FC_compareGoBp = simplify(FC_compareGoBp)

#
summary(simple_LA_compareGoBp)[summary(simple_LA_compareGoBp)$ID %in% intersect(summary(simple_LA_compareGoBp)$ID, summary(simple_FC_compareGoBp)$ID),'Description'] #simplify then intersect
 
#save shared GoBp
shared_GoBp = summary(simple_LA_compareGoBp)[summary(simple_LA_compareGoBp)$ID %in% intersect(summary(LA_compareGoBp)$ID, summary(FC_compareGoBp)$ID),] # intersect then simplify
shared_GoBp = shared_GoBp[,c(c(2,3,1),4:10)]
colnames(shared_GoBp)[3:10] <- paste0("LA_", colnames(shared_GoBp)[3:10])

FC_compareGoBp_merge = summary(FC_compareGoBp)
FC_compareGoBp_merge = FC_compareGoBp_merge[match(shared_GoBp$ID,FC_compareGoBp_merge$ID),]
colnames(FC_compareGoBp_merge) = paste0("FC_", colnames(FC_compareGoBp_merge) )
shared_GoBp = cbind(shared_GoBp, FC_compareGoBp_merge[,c(1,(4:10) )])
write.csv(shared_GoBp, file = 'csvs/GO_SharedGoBp_Injury_SST-RMA_R_1e-2.csv',row.names=FALSE)

