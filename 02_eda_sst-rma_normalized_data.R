## 02 - Exploratory data analysis and plotting for normalized data
# Normalization done with SST-RMA using all 16 samples (both FC and LA; all three conditions together).
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/normalized_data_n16_sst-rma.rda')
library(jaffelab)
library(ggplot2)
theme_set(theme_bw(base_size=40) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))
## PCA
pca = prcomp(t(sst_rma_signal),scale.=TRUE)
varPCs=getPcaVars(pca)
exprsPCs=pca$x
pd=cbind(pd,exprsPCs[,1:10])
pd$PC_VarExplained = varPCs

# Plotting
PC1_2_region = ggplot(data=pd, aes(x=PC1,y=PC2,col=brain_region ) ) + 
  geom_point(size=5) +	 
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x=paste0("PC1"," (", as.character(signif(pd$PC_VarExplained[1],3)),"%)" ),
       y=paste0("PC2"," (", as.character(signif(pd$PC_VarExplained[2],3)),"%)" )) +
  theme(legend.position = c(.9,.15), 
        legend.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Dark2")

PC1_2_group = ggplot(data=pd, aes(x=PC1,y=PC2,col=Group ) ) + 
  geom_point(size=5) +	 
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") + 
  labs(x=paste0("PC1"," (", as.character(signif(pd$PC_VarExplained[1],3)),"%)" ),
       y=paste0("PC2"," (", as.character(signif(pd$PC_VarExplained[2],3)),"%)" )) +
  theme(legend.position = c(.85,.15), 
        legend.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines') )+ scale_colour_brewer(palette="Set1")

pdf('plots/PCA_SST-RMA_Normalization.pdf', height=12,width=12)
PC1_2_region
PC1_2_group
dev.off()

tiff('plots/PCA_SST-RMA_Normalization_S1A.tiff', height=12,width=12, units = "in", res = 600)
PC1_2_region
dev.off()
tiff('plots/PCA_SST-RMA_Normalization_S1B.tiff', height=12,width=12, units = "in", res = 600)
PC1_2_group
dev.off()



# PCA statistics
summary(lm(pd$PC1~ (as.factor(pd$Group))+pd$brain_region ))
summary(aov(pd$PC1~pd$Group))
t.test(pd$PC1~pd$brain_region)

summary(lm(pd$PC2~ (as.factor(pd$Group))+pd$brain_region ))
summary(aov(pd$PC2~pd$Group))
t.test(pd$PC2~pd$brain_region)
