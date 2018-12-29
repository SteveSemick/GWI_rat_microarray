### Open .chp files from TAC
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')

### Compare results from RMA normalization with subsequent testing 

# Load R Stats
load('rda/mergedStats_RMA_R.rda',verbose=T) 

# Load TAC STats
load('TAC/RMA_SingleRegion/MergedStats.rda',verbose=T) # res_stats

res_stats_P = sapply(res_stats, function(x) x[match(rownames(mergedStats) , x[,'ID']),'P-val'])
rownames(res_stats_P) = rownames(mergedStats)
colnames(res_stats_P) <- paste0(colnames(res_stats_P), "_P")

##
cor.test(res_stats_P[,'FC_AvB_P'], mergedStats[,'FC_AvB_P.Value'])
cor.test(res_stats_P[,'FC_AvC_P'], mergedStats[,'FC_AvC_P.Value'])
cor.test(res_stats_P[,'FC_BvC_P'], mergedStats[,'FC_BvC_P.Value'])
#
cor.test(res_stats_P[,'LA_AvB_P'], mergedStats[,'LA_AvB_P.Value'])
cor.test(res_stats_P[,'LA_AvC_P'], mergedStats[,'LA_AvC_P.Value'])
cor.test(res_stats_P[,'LA_BvC_P'], mergedStats[,'LA_BvC_P.Value'])


## FC ##
pdf('plots/Comparison_RMA_Model_Pvalues_TAC_vs_R.pdf', height=12, width=12)
par(mar=c(5,6,2,2))
plot(-log10(res_stats_P[,'FC_AvB_P']), -log10(mergedStats[,'FC_AvB_P.Value']), pch = 21, bg="grey",
     xlab="-log10(P) [From TAC]", ylab="-log10(P) [From R]", main ="FC: AvB" ,
     xlim=c(0,5), ylim = c(0,5),
     cex.axis=2,cex.lab=2, cex.main=2)
lines(x = c(0,5), y = c(0,5),col="red")
##
plot(-log10(res_stats_P[,'FC_AvC_P']), -log10(mergedStats[,'FC_AvC_P.Value']), pch = 21, bg="grey",
     xlab="-log10(P) [From TAC]", ylab="-log10(P) [From R]", main ="FC: AvC" ,
     xlim=c(0,5), ylim = c(0,5),
     cex.axis=2,cex.lab=2, cex.main=2)
lines(x = c(0,5), y = c(0,5),col="red")
##
plot(-log10(res_stats_P[,'FC_BvC_P']), -log10(mergedStats[,'FC_BvC_P.Value']), pch = 21, bg="grey",
     xlab="-log10(P) [From TAC]", ylab="-log10(P) [From R]", main ="FC: BvC" ,
     xlim=c(0,5), ylim = c(0,5),
     cex.axis=2,cex.lab=2, cex.main=2)
lines(x = c(0,5), y = c(0,5),col="red")

## LA ##
plot(-log10(res_stats_P[,'LA_AvB_P']), -log10(mergedStats[,'LA_AvB_P.Value']), pch = 21, bg="grey",
     xlab="-log10(P) [From TAC]", ylab="-log10(P) [From R]", main ="LA: AvB" ,
     xlim=c(0,5), ylim = c(0,5),
     cex.axis=2,cex.lab=2, cex.main=2)
lines(x = c(0,5), y = c(0,5),col="red")
##
plot(-log10(res_stats_P[,'LA_AvC_P']), -log10(mergedStats[,'LA_AvC_P.Value']), pch = 21, bg="grey",
     xlab="-log10(P) [From TAC]", ylab="-log10(P) [From R]", main ="LA: AvC" ,
     xlim=c(0,5), ylim = c(0,5),
     cex.axis=2,cex.lab=2, cex.main=2)
lines(x = c(0,5), y = c(0,5),col="red")
##
plot(-log10(res_stats_P[,'LA_BvC_P']), -log10(mergedStats[,'LA_BvC_P.Value']), pch = 21, bg="grey",
     xlab="-log10(P) [From TAC]", ylab="-log10(P) [From R]", main ="LA: BvC" ,
     xlim=c(0,5), ylim = c(0,5),
     cex.axis=2,cex.lab=2, cex.main=2)
lines(x = c(0,5), y = c(0,5),col="red")
dev.off()
#res_stats[[1]][match(rownames(mergedStats), res_stats[[1]][,'ID'] ),'ID'] == rownames(mergedStats)
