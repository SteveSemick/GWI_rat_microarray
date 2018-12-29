setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
## Required packages
require(openxlsx)
require(data.table)

################ Both-region models ################
## Import simple model results from 
RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/RMA_BothRegion',pattern = '.txt',full.names = T)

res = lapply(RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(RMA_Res))


### Extract Normalized Signal
rma_signal = res[[1]]
rownames(rma_signal) <- rma_signal$ID
rma_signal = rma_signal[, grep("Signal", colnames(rma_signal)) ]

### Extract Statistics from each model
res_stats = lapply(res, function(x) x[,1:21] )
res_stats = lapply(res_stats, function(x) x[,-c(11,12)] )
res_stats = lapply(res_stats, function(x) x[order(x[,'P-val']),] ) 

###
res_FDR05 = lapply(res_stats, function(x) x[x[,'FDR P-val']<0.05,])

### Export results

# Write rda
save(res_stats, file='TAC/RMA_BothRegion/MergedStats.rda')
save(rma_signal, file='TAC/RMA_BothRegion/Norm_signal.rda')

################ Single-region models ################
## Import simple model results from 
RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/RMA_SingleRegion',pattern = '.txt',full.names = T)

res = lapply(RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(RMA_Res))


### Extract Normalized Signal
rma_signal = res[[1]]
rownames(rma_signal) <- rma_signal$ID
rma_signal = rma_signal[, grep("Signal", colnames(rma_signal)) ]

### Extract Statistics from each model
res_stats = lapply(res, function(x) x[,1:21] )
res_stats = lapply(res_stats, function(x) x[,-c(11,12)] )
res_stats = lapply(res_stats, function(x) x[order(x[,'P-val']),] ) 

###
res_FDR05 = lapply(res_stats, function(x) x[x[,'FDR P-val']<0.05,])

### Export results

# Write rda
save(res_stats, file='TAC/RMA_SingleRegion/MergedStats.rda')
save(rma_signal, file='TAC/RMA_SingleRegion/Norm_signal.rda')

