### Summarise results from TAC SST-RMA Full model

load('TAC/SST-RMA_BothRegions/MergedStats.rda')
names(res_stats) <- gsub("_v_", "v", names(res_stats) )

merged_SstRma = res_stats[[1]]
colnames(merged_SstRma)[(2:11)] <- paste0( colnames(merged_SstRma)[(2:11)], "_", names(res_stats)[1] )

tmp = lapply(res_stats[2:6], function(x) {x[match(merged_SstRma$ID,x[,'ID']),2:11]} )
newNames = lapply(names(tmp), function(x) paste0(colnames(merged_SstRma)[2:11], "_",x ) )
colnames(tmp[[1]]) <- newNames[[1]]
colnames(tmp[[2]]) <- newNames[[2]]
colnames(tmp[[3]]) <- newNames[[3]]
colnames(tmp[[4]]) <- newNames[[4]]
colnames(tmp[[5]]) <- newNames[[5]]

#tmp2  = mapply( function(x,y) setNames(object=x, y), tmp, newNames)
merged_SstRma = cbind(merged_SstRma, do.call("cbind", tmp))
