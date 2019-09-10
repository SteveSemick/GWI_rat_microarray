library(ggplot2)
theme_set(theme_bw(base_size=20) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"))

## Other post-hoc tests
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/injury_DE_res.rda', verbose=T)

#Altered inflammatory activity associated with reduced hippocampal volume and more severe posttraumatic stress symptoms in Gulf War veterans.
merged_SstRma_disease[grep('Tnfrsf1b',merged_SstRma_disease$Symbol),]



dat = merged_SstRma_disease[grep('cholinergic receptor',as.character(merged_SstRma_disease$gene_description)),]

dat$postHoc_NomSig = dat[,'LA_CTRL_v_PYR_LPS-P.Value']<0.05

dat$Symbol = factor(dat$Symbol, levels = dat$Symbol[order(dat[,'LA_CTRL_v_PYR_LPS-logFC'])])

ACh_receptor_plot = ggplot(dat,aes(x=Symbol,y=`LA_CTRL_v_PYR_LPS-logFC`,fill=postHoc_NomSig)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
  labs(x='ACh Receptor', 
       y='log2 Fold Change', 
       title='ACh Receptor Expression \nAfter Injury in Lateral Amygdala') + 
geom_text(aes(x=Symbol, 
              y= -0.015*sign(`LA_CTRL_v_PYR_LPS-logFC`), 
              label=Symbol, 
              colour=postHoc_NomSig,
              vjust=ifelse(sign(`LA_CTRL_v_PYR_LPS-logFC`)<0, 0.5, 0.5),
          hjust=ifelse(sign(`LA_CTRL_v_PYR_LPS-logFC`)<0, 0, 1),
          angle=ifelse(sign(`LA_CTRL_v_PYR_LPS-logFC`)<0, 60, 60) 
          ),size=5) + scale_fill_manual(values=c("black","#ab1323" ) ) +
scale_colour_manual(values=c("black","#ab1323" ) )
ggsave(ACh_receptor_plot, file = "plots/Acetlycholine Receptor Post-Injury in LA.tiff",
       height=8,width=8,units="in", dpi=600)