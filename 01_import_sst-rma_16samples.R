## 01 - Import normalized data from TAC
# Normalization done with SST-RMA using all 16 samples (both FC and LA; all three conditions together).
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
## Import full  model results from 
SST_RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/SST-RMA_BothRegions', pattern = '.txt', full.names = T)

res = lapply(SST_RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(SST_RMA_Res))

### Extract Normalized Signal
sst_rma_signal = res[[1]]
rownames(sst_rma_signal) <- sst_rma_signal$ID
sst_rma_signal = sst_rma_signal[, grep("Signal", colnames(sst_rma_signal)) ]

## Create phenotype table
pd = data.frame(fileNames= colnames(sst_rma_signal) )
pd$sampleNames <- gsub("\\..*","", pd$fileNames)
  
pd$brain_region = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",1)

pd$sampleNum = jaffelab::ss(jaffelab::ss(pd$sampleNames,"_",2),"-",2)
pd$sampleNum = as.numeric(gsub("\\..*","",pd$sampleNum))

pd$Group = ifelse(pd$sampleNum %in% c(1,2,3), "CTRL", ifelse(pd$sampleNum %in% c(4,5,6), "PYR+LPS", "PYR+LPS+RGZ") ) # assigning
#CTRL - Control
#LPS - lipopolysacharide
#PYR - pyridostigmine
#RGZ - roglitizone
pd=pd[,-1]

#Fix column names for normalized signal data
colnames(sst_rma_signal) <- gsub("\\..*","", colnames(sst_rma_signal))
save(pd, sst_rma_signal, file='rda/normalized_data_n16_sst-rma.rda')
