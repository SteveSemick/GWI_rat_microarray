## RMA Normalization and QC
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
library(affycoretools)
theme_set(theme_bw(base_size=20) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))

load('rda/twoBatch_rawData.rda')

###=========== Normalization =====###
normDat = rma(rawDat_batchBoth) # rma normalization
normDat = getMainProbes(normDat) # drop control probes 

normExprs = exprs(normDat)
pd = pdM = pData(normDat)

save(pd,normExprs, file='rda/twoBatch_rmaNorm_Data.rda')
###=========== Quality Control =====###

myColors <- ifelse(pdM$Batch=="Batch1" , rgb(0.1,0.1,0.7,0.5), rgb(0.8,0.1,0.3,0.6))

pdf('plots/QC_Batch_normExpression_Boxplot.pdf')
boxplot(normExprs, ylim=c(0, 14), col = myColors,ylab='log2 intensity',xlab='array',xaxt ='n',outline=F)
dev.off()

rma_medians <- rowMedians((normExprs))
pdf('plots/histogram_of_postRMA_median_intensity.pdf')
hist(rma_medians, 100, col = "cornsilk", freq = FALSE, main = "Histogram of median intensities", border = "antiquewhite4", xlab = "Median intensities") 
dev.off()

### PCA

pca = prcomp(t(normExprs),scale.=F)
varPCs=getPcaVars(pca)
exprsPCs=pca$x
pd=cbind(pdM,exprsPCs[,1:5])
pd$PC_VarExplained = varPCs


# Plotting
t.test(PC1~Batch,data=pd)
boxplot(PC1~Batch,data=pd)
boxplot(PC1~brain_region,data=pd)

PC1_2_batch = ggplot(data=pd, aes(x=PC1,y=PC2,col=Batch,shape=brain_region ) ) + 
  geom_point(size=3) +	 
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  scale_shape_manual(values = c(15,17)) +
  labs(x=paste0("PC1"," (", as.character(signif(pd$PC_VarExplained[1],3)),"%)" ),
       y=paste0("PC2"," (", as.character(signif(pd$PC_VarExplained[2],3)),"%)" )) + theme(legend.position='right')
pdf('plots/PCA_rmaNorm_Data_PC1_2_Batch_X_Region.pdf')
PC1_2_batch
dev.off()
#theme(legend.position = c(.9,.15), 
#     legend.background = element_rect(colour = "black"),
#    legend.title=element_blank(),
#   legend.key = element_rect(size = 5),
#  legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2")
PC1_2_batch
##

t.test(pd$PC2,pd$colMed)
boxplot(PC2~Batch,data=pd)
##
pd$Group=factor(pd$Group,levels=c('Naive','CTRL','PYR_LPS','PYR_LPS_RGZ'))
t.test(PC3~brain_region,data=pd)
boxplot(PC3~Group,data=pd)

boxplot(PC2~Group,data=pd)

boxplot(PC4~Group,data=pd)
boxplot(PC4~brain_region,data=pd)
boxplot(PC4~Batch,data=pd)
boxplot(PC5~Group,data=pd)

boxplot(PC3~Group,data=pd)
summary(aov(PC1~Group,data=pd))
summary(aov(PC2~Group,data=pd))
summary(aov(PC3~Group,data=pd))
summary(aov(PC4~Group,data=pd))
summary(aov(PC5~Group,data=pd))
