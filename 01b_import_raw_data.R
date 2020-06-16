#Script to import and analyze PTSD Rat Clarion S arrays
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
##
library(limma)
library(GEOquery)
library(clariomsrattranscriptcluster.db)
## Import .cel files
library(oligo)
library(ggplot2)
library(jaffelab)
library(plyr)
library(tidyr)
theme_set(theme_bw(base_size=40) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))
## Import cels
cel_path1 = "C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/RatClariomS_16"
celFiles_batch1 = list.files(cel_path1,full.names=TRUE,pattern =".CEL")

cel_path2="C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/Woo_RatClariomS_18_191011"
celFiles_batch2 = list.files(cel_path2,full.names=TRUE,pattern =".CEL")

#rawDat_batch1 = read.celfiles(celFiles_batch1)


rawDat_batchBoth = read.celfiles(celFiles_batch1,celFiles_batch2)

## Create phenotype table
pd = phenoData(rawDat_batchBoth)
pd = as(pd, "data.frame")
pd$sampleNames = rownames(pd)

pd$Batch = ifelse(grepl("Woo-RatClariomS-18",pd$sampleNames),"Batch2","Batch1")

## Process batch 1 phenotype data
pd1 = pd[pd$Batch=="Batch1",]

pd1$brain_region = jaffelab::ss(jaffelab::ss(pd1$sampleNames,"_",2),"-",1)
##
pd1$sampleNum = jaffelab::ss(jaffelab::ss(pd1$sampleNames,"_",2),"-",2)
pd1$sampleNum = as.numeric(gsub("\\..*","",pd1$sampleNum))
##
pd1$Group = ifelse(pd1$sampleNum %in% c(1,2,3), "CTRL", ifelse(pd1$sampleNum %in% c(4,5,6), "PYR_LPS", "PYR_LPS_RGZ") ) # assigning
##
pd1$RatNum =NA
pd1$Batch="Batch1"
pd1=pd1[,c('sampleNames','sampleNum','RatNum','Batch','brain_region','Group')]

## Process batch 2 phenotype data
## Create phenotype table
pd2 = pd[pd$Batch=="Batch2",]

pd2$sampleNum = jaffelab::ss(jaffelab::ss(pd2$sampleNames,"_",2),"-",1)
pd2$sampleNum = as.numeric(gsub("\\..*","",pd2$sampleNum))

## Getting the rat number
pd2$RatNum = NA
pd2[pd2$sampleNum%in%(1:9),'RatNum'] = pd2[pd2$sampleNum%in%(1:9),'sampleNum']+51
pd2[pd2$sampleNum%in%(10:18),'RatNum'] = pd2[pd2$sampleNum%in%(10:18),'sampleNum']+42

## Adding brain region
pd2$brain_region = NA
pd2$brain_region[pd2$sampleNum%in%(1:9)] = "FC"
pd2$brain_region[pd2$sampleNum%in%(10:18)] = "LA"

## Batch2 group
pd2$Group = NA
pd2$Group[pd2$RatNum%in%c(52,53,54)] = "PYR_LPS"
pd2$Group[pd2$RatNum%in%c(55,56,57)] = "CTRL"
pd2$Group[pd2$RatNum%in%c(58,59,60)] = "Naive"

#
pd2$Batch = "Batch2"
pd2=pd2[,c('sampleNames','sampleNum','RatNum','Batch','brain_region','Group')]

## Merge
pdM=rbind(pd1,pd2)
pData(rawDat_batchBoth) = pdM
pdM= pdM[order(pdM$Batch,pdM$sampleNum),]
write.csv(pdM,file='csvs/pd_n34.csv',row.names=FALSE)
save(rawDat_batchBoth,file='rda/twoBatch_rawData.rda')
