### Summarise results from TAC SST-RMA Full model

load('TAC/SST-RMA_BothRegions/MergedStats.rda')
names(res_stats) <- gsub("_v_", "v", names(res_stats) )

merged_SstRma = res_stats[[1]]
colnames(merged_SstRma)[(2:11)] <- paste0( names(res_stats)[1], ".", colnames(merged_SstRma)[(2:11)]  )

tmp = lapply(res_stats[2:6], function(x) {x[match(merged_SstRma$ID,x[,'ID']),2:11]} )
#newNames = lapply(names(tmp), function(x) paste0(colnames(merged_SstRma)[2:11], "_",x ) )
#colnames(tmp[[1]]) <- newNames[[1]]
#colnames(tmp[[2]]) <- newNames[[2]]
#colnames(tmp[[3]]) <- newNames[[3]]
#colnames(tmp[[4]]) <- newNames[[4]]
#colnames(tmp[[5]]) <- newNames[[5]]

#tmp2  = mapply( function(x,y) setNames(object=x, y), tmp, newNames)
merged_SstRma = cbind(merged_SstRma, do.call("cbind", tmp))
save(merged_SstRma, file='rda/merged_sstRma_results.rda')

FDR_Col = grep("FDR",colnames(merged_SstRma),value=T)

#### 
SigNums = data.frame(colSums(merged_SstRma[,FDR_Col]<0.10))
SigNums$Region = gsub("_.*","", rownames(SigNums))

################ Making venn diagram at Probe level #########
SigIDs = sapply(FDR_Col, function(x) merged_SstRma[merged_SstRma[,x]<0.10,'ID'])

library(VennDiagram)
## sig FC
FC_Venn_ID = venn.diagram(x = SigIDs[grep("FC", names(SigIDs),value=T)],
                                 category.names = c('A v. B','A v. C', 'B v. C'),
                                 filename = NULL,
                                 fill = c('red', 'blue','green'),
                                 #cat.just=list(c(0.9,1.5) , c(-0.8,5) ), 
                                 cex=5, cat.cex=2)

pdf('plots/FC_Venn_ID_FDR10.pdf')
grid.draw(FC_Venn_ID)
dev.off()

## sig LA
LA_Venn_ID = venn.diagram(x = SigIDs[grep("LA", names(SigIDs),value=T)],
                          category.names = c('A v. B','A v. C', 'B v. C'),
                          filename = NULL,
                          fill = c('red', 'blue','green'),
                          #cat.just=list(c(0.9,1.5) , c(-0.8,5) ), 
                          cex=5, cat.cex=2)

## sig AvB
AvB_Venn_ID = venn.diagram(x = SigIDs[grep("AvB", names(SigIDs),value=T)],
                          category.names = c('FC','LA'),
                          filename = NULL,
                          fill = c('red', 'blue'),
                          main = 'A vs. B',
                          main.cex = 5,
                          #cat.just=list(c(0.9,1.5) , c(-0.8,5) ), 
                          cex=5, cat.cex=2)

## sig AvC
AvC_Venn_ID = venn.diagram(x = SigIDs[grep("AvC", names(SigIDs),value=T)],
                          category.names = c('FC','LA'),
                          filename = NULL,
                          fill = c('red', 'blue'),
                          main = 'A vs. C',
                          main.cex = 5,
                          #cat.just=list(c(0.9,1.5) , c(-0.8,5) ), 
                          cex=5, cat.cex=2)
## sig BvC
BvC_Venn_ID = venn.diagram(x = SigIDs[grep("BvC", names(SigIDs),value=T)],
                          category.names = c('FC','LA'),
                          filename = NULL,
                          fill = c('red', 'blue'),
                          main = 'B vs. C',
                          main.cex = 5,
                          #cat.just=list(c(0.9,1.5) , c(-0.8,5) ), 
                          cex=5, cat.cex=2)
pdf('plots/Comparison_Overlap_BetweenRegions_Venn_ID_FDR10.pdf')
grid.draw(AvB_Venn_ID)
grid.newpage()
grid.draw(AvC_Venn_ID)
grid.newpage()
grid.draw(BvC_Venn_ID)
dev.off()
