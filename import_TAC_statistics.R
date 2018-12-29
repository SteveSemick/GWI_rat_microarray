setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
## Required packages
require(openxlsx)
require(data.table)

## Import simple model results from 
SST_RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/SST-RMA_BothRegions',pattern = '.txt',full.names = T)

res = lapply(SST_RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(SST_RMA_Res))


### Extract Normalized Signal
sst_rma_signal = res[[1]]
rownames(sst_rma_signal) <- sst_rma_signal$ID
sst_rma_signal = sst_rma_signal[, grep("Signal", colnames(sst_rma_signal)) ]

### Extract Statistics from each model
res_stats = lapply(res, function(x) x[,1:21] )
res_stats = lapply(res_stats, function(x) x[,-c(11,12)] )
res_stats = lapply(res_stats, function(x) x[order(x[,'P-val']),] ) 

###
res_FDR05 = lapply(res_stats, function(x) x[x[,'FDR P-val']<0.05,])

### Export results
# Write xlsx file
openxlsx::write.xlsx(res_stats, file='TAC/SST-RMA_BothRegions/MergedStats.xlsx' )
openxlsx::write.xlsx(res_FDR05, file='TAC/SST-RMA_BothRegions/MergedStats_FDR05.xlsx' )

# Write rda
save(res_stats, file='TAC/SST-RMA_BothRegions/MergedStats.rda')
save(sst_rma_signal, file='TAC/SST-RMA_BothRegions/Norm_signal.rda')

