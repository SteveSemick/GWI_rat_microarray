## Plotting fold-change versus normalized transcript values
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
library(matrixStats)
## Load stats
load('rda/merged_sstRma_results.rda',verbose=T)

## Extract normalized values

RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/SST-RMA_BothRegions',pattern = '.txt',full.names = T)

res = lapply(RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(RMA_Res))

### Extract Normalized Signal
rma_signal = res[[1]]
rownames(rma_signal) <- rma_signal$ID
rma_signal = rma_signal[, grep("Signal", colnames(rma_signal)) ]

pd=data.frame(sampleNames=colnames(rma_signal),stringsAsFactors=FALSE)
pd$brain_region = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",1)

pd$sampleNum = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",2)
pd$sampleNum = as.numeric(gsub("\\..*","",pd$sampleNum))

pd$Group = ifelse(pd$sampleNum %in% c(1,2,3), "GroupA", ifelse(pd$sampleNum %in% c(4,5,6), "GroupB", "GroupC") ) # assigning


### Plotting
rownames(merged_SstRma) <- merged_SstRma$ID
hist(rowMeans(rma_signal[ , pd[pd$Group%in%c('GroupA','GroupB')&pd$brain_region=="FC",'sampleNames'] ] ) )

plot(log2(abs(merged_SstRma[rownames(rma_signal),'FC_AvB.Fold Change'] ) ), rowMeans(rma_signal[ , pd[pd$Group%in%c('GroupA','GroupB')&pd$brain_region=="FC",'sampleNames'] ] ) )  
plot(rowSds(as.matrix(rma_signal[ , pd[pd$Group%in%c('GroupA','GroupB')&pd$brain_region=="FC",'sampleNames'] ] ) ), rowMeans(rma_signal[ , pd[pd$Group%in%c('GroupA','GroupB')&pd$brain_region=="FC",'sampleNames'] ] ) )  

## All genes