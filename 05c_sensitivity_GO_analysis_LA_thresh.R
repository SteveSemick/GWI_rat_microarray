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

################## MORE STRINGENT P-VALUE THRESHOLD (p<0.001 ##################

## up and down stats - LA
sigStats = merged_SstRma_disease[which(merged_SstRma_disease[,'LA_CTRL_v_PYR_LPS-P.Value'] < 1e-3), c("LA_CTRL_v_PYR_LPS-logFC", "LA_CTRL_v_PYR_LPS-adj.P.Val", "EntrezID")]
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
#LA_compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Rn.eg.db, qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE) #n.s.

LA_compareGo = lapply(list(compareGoMf=LA_compareGoMf, compareGoBp=LA_compareGoBp),simplify)

## Export results
save(LA_compareGo,LA_compareKegg, LA_compareGoMf,LA_compareGoBp, file='rda/LA_GO_Results_sensitivity_1e-3_split.rda')
openxlsx::write.xlsx(lapply(c(KEGG=LA_compareKegg, LA_compareGo), summary), file = 'csvs/GO_LA_Injury_SST-RMA_R_1e-3_split.xlsx')

################## LESS STRINGENT P-VALUE THRESHOLD (p<0.05 ##################

## up and down stats - LA
sigStats = merged_SstRma_disease[which(merged_SstRma_disease[,'LA_CTRL_v_PYR_LPS-P.Value'] < 0.05), c("LA_CTRL_v_PYR_LPS-logFC", "LA_CTRL_v_PYR_LPS-adj.P.Val", "EntrezID")]
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

LA_compareGo = lapply(list(compareGoMf=LA_compareGoMf, compareGoBp=LA_compareGoBp),simplify)

## Export results
save(LA_compareGo,LA_compareKegg, LA_compareGoMf,LA_compareGoBp, file='rda/LA_GO_Results_sensitivity_05_split.rda')
openxlsx::write.xlsx(lapply(c(KEGG=LA_compareKegg, LA_compareGo), summary), file = 'csvs/GO_LA_Injury_SST-RMA_R_05_split.xlsx')

## Finding the overlap between the different thresholds

## Check for overlap the different LA conditions
load('rda/LA_GO_Results.rda',verbose=T)
GoKegg_01 = LA_compareKegg
GoMf_01 =LA_compareGoMf
GoBp_01 = LA_compareGoBp

load('rda/LA_GO_Results_sensitivity_1e-3_split.rda',verbose=T) # More stringent threshold
GoKegg_001 = LA_compareKegg
GoMf_001 =LA_compareGoMf
GoBp_001 = LA_compareGoBp

load('rda/LA_GO_Results_sensitivity_05_split.rda',verbose=T) #Less stringent threshold
GoKegg_05 = LA_compareKegg
GoMf_05 =LA_compareGoMf
GoBp_05 = LA_compareGoBp

Reduce(intersect, list(summary(GoKegg_01)$ID, summary(GoKegg_001)$ID, summary(GoKegg_05)$ID) )
Reduce(intersect, list(summary(GoMf_01)$ID, summary(GoMf_001)$ID, summary(GoMf_05)$ID) )
Reduce(intersect, list(summary(GoBp_01)$ID, summary(GoBp_001)$ID, summary(GoBp_05)$ID) )

simple_01_compareGoMf = simplify(GoMf_01)
simple_01_compareGoBp = simplify(GoBp_01)

#
summary(simple_LA_compareGoBp)[summary(simple_LA_compareGoBp)$ID %in% intersect(summary(simple_LA_compareGoBp)$ID, summary(simple_FC_compareGoBp)$ID),'Description'] #simplify then intersect

#save shared GoMf
shared_GoMf = summary(simple_01_compareGoMf)
shared_GoMf$Sig_Sensitivity = summary(simple_01_compareGoMf)$ID %in% Reduce(intersect, list(summary(GoMf_01)$ID, summary(GoMf_001)$ID, summary(GoMf_05)$ID) ) # intersect then simplify
shared_GoMf = shared_GoMf[,c(c(2,3,1),4:11)]

#save shared GoBp
shared_GoBp = summary(simple_01_compareGoBp)
shared_GoBp$Sig_Sensitivity = summary(simple_01_compareGoBp)$ID %in% Reduce(intersect, list(summary(GoBp_01)$ID, summary(GoBp_001)$ID, summary(GoBp_05)$ID) ) # intersect then simplify
shared_GoBp = shared_GoBp[,c(c(2,3,1),4:11)]

FC_compareGoBp_merge = summary(FC_compareGoBp)
FC_compareGoBp_merge = FC_compareGoBp_merge[match(shared_GoBp$ID,FC_compareGoBp_merge$ID),]
colnames(FC_compareGoBp_merge) = paste0("FC_", colnames(FC_compareGoBp_merge) )
shared_GoBp = cbind(shared_GoBp, FC_compareGoBp_merge[,c(1,(4:10) )])
write.csv(shared_GoBp, file = 'csvs/GO_SharedGoBp_Injury_SST-RMA_R_1e-2.csv',row.names=FALSE)


### Visualizations
cnetplot(simple_01_compareGoBp, foldChange=geneList)
enrichMap(shared_GoBp, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)


### WikiPathways Analysis
library(magrittr)
library(clusterProfiler)

## up and down stats - LA
sigStats = merged_SstRma_disease[which(merged_SstRma_disease[,'LA_CTRL_v_PYR_LPS-P.Value'] < 0.05), c("LA_CTRL_v_PYR_LPS-logFC", "LA_CTRL_v_PYR_LPS-adj.P.Val", "EntrezID")]
sigStats$Sign = sign(sigStats$`LA_CTRL_v_PYR_LPS-logFC`)
sigStats = sigStats[!is.na(sigStats$EntrezID),]
sigStats = sigStats[!duplicated(sigStats[,c("EntrezID")]),]
sigStats =sigStats[-grep(";",sigStats$EntrezID),]
sigStats$EntrezID = as.character(sigStats$EntrezID)
## by sign
gList = split(sigStats$EntrezID, sigStats$Sign)
lengths(gList)
gList = lapply(gList,as.character)
###

download.file('http://data.wikipathways.org/current/gmt/wikipathways-20190810-gmt-Rattus_norvegicus.gmt', destfile = 'C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/wikipathway_20190810.gmt')
wpgmtfile <- 'C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/wikipathway_20190810.gmt'
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewpDown <- enricher(gList[['-1']], TERM2GENE = wpid2gene, TERM2NAME = wpid2name, universe=univ,pvalueCutoff = 1, qvalueCutoff = 1)
ewpUp <- enricher(gList[['1']], TERM2GENE = wpid2gene, TERM2NAME = wpid2name, universe=univ, pvalueCutoff=1)

sigUnsplit = merged_SstRma_disease[which(merged_SstRma_disease[,'LA_CTRL_v_PYR_LPS-P.Value'] < 0.001), c( "EntrezID")]
sigUnsplit = sigUnsplit[!is.na(sigUnsplit)]
sigUnsplit = sigUnsplit[!duplicated(sigUnsplit)]
sigUnsplit =sigUnsplit[-grep(";",sigUnsplit)]
sigUnsplit= as.character(sigUnsplit)
ewpUnsplit =enricher(sigUnsplit, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, universe=univ, pvalueCutoff=0.2,qvalue=0.2) 
emapplot(simple_01_compareGoBp)

