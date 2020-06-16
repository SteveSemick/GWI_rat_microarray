## QC on raw data #######
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
theme_set(theme_bw(base_size=20) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))
load('rda/twoBatch_rawData.rda',verbose=T)

### Boxplots
rawExprs = exprs(rawDat_batchBoth)
pdM=pData(rawDat_batchBoth)
# Prepare a vector of colors with specific color for Nairobi and Eskimo
myColors <- ifelse(pdM$Batch=="Batch1" , rgb(0.1,0.1,0.7,0.5), rgb(0.8,0.1,0.3,0.6))

pdf('plots/QC_Batch_RawExpression_Boxplot.pdf')
boxplot(log2(rawExprs), ylim=c(3, 13), col = myColors,ylab='log2 intensity',xlab='array',xaxt ='n',outline=F)
dev.off()

### PCA
pca = prcomp(t(log2(rawExprs)),scale.=FALSE)
varPCs=getPcaVars(pca)
exprsPCs=pca$x
pd=cbind(pdM,exprsPCs[,1:5])
pd$PC_VarExplained = varPCs


# Plotting
t.test(PC1~Batch,data=pd)
boxplot(PC1~Batch,data=pd)

PC1_2_batch = ggplot(data=pd, aes(x=PC1,y=PC2,col=Batch ) ) + 
  geom_point(size=3) +	 
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x=paste0("PC1"," (", as.character(signif(pd$PC_VarExplained[1],3)),"%)" ),
       y=paste0("PC2"," (", as.character(signif(pd$PC_VarExplained[2],3)),"%)" )) +
  theme(legend.position = 'right')
pdf('plots/PCA_rawData_PC1_2_Batch.pdf')
PC1_2_batch
dev.off()

##
pd$colMed=matrixStats::colMedians(as.matrix(log2(rawExprs)))

t.test(pd$PC2,pd$colMed)
boxplot(PC2~Batch,data=pd)
##
t.test(PC3~brain_region,data=pd)
boxplot(PC3~brain_region,data=pd)