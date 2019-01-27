### Check post-hoc hypothesis about PPAR-gamma related genes
setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')

load('rda/merged_sstRma_results.rda',verbose=T)
library(dplyr)
library(ggplot2)
library('ggrepel')
theme_set(theme_bw(base_size=24) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"				 )) 


### Calculating LFC
merged_SstRma$LA_AvB_log2FC = merged_SstRma[,'LA_AvB.LA: A Avg (log2)'] - merged_SstRma[,'LA_AvB.LA: B Avg (log2)']
merged_SstRma$FC_AvB_log2FC = merged_SstRma[,'FC_AvB.FC: A Avg (log2)'] - merged_SstRma[,'FC_AvB.FC: B Avg (log2)']

merged_SstRma$LA_BvC_log2FC = merged_SstRma[,'LA_BvC.LA: B Avg (log2)'] - merged_SstRma[,'LA_BvC.LA: C Avg (log2)']
merged_SstRma$FC_BvC_log2FC = merged_SstRma[,'FC_BvC.FC: B Avg (log2)'] - merged_SstRma[,'FC_BvC.FC: C Avg (log2)']


#### Genes from Villapol 2018
human_PPARG_gene = c("PPARG","RXRA", "RXRG", "RXRB", "CD36", "JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","FOSL1P1","STAT1", "STAT3", "STAT6", "NFKB1","NFKB2","REL","RELA","RELB", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "MMP1", "MMP3", "MMP9", "IL10", "ARG1", "NOS2","PTGS2","IL6", "TNF","SOD1","SOD2","SOD3")  


### Map human genes found in literature to rat orthologs
library(biomaRt)
ensembl <- useMart("ensembl")
human = useDataset("hsapiens_gene_ensembl", mart=ensembl)
rat = useDataset("rnorvegicus_gene_ensembl", mart=ensembl)


humanAttributes = c("ensembl_gene_id","entrezgene","hgnc_symbol")
ratAttributes = c("ensembl_gene_id","entrezgene","rgd_symbol","description","chromosome_name","start_position","end_position","strand","gene_biotype","ensembl_transcript_id")
### query biomart
#ratAttributes = searchAttributes(mart = human, pattern = "rnor")
#ratMapped_PPARG_Genes = getBM(attributes = c('hgnc_symbol',c(ratAttributes$name) ),
#      filters = 'hgnc_symbol', 
#      values = human_PPARG_gene, 
#      mart = human)

ratMapped_PPARG_Genes = getLDS(attributes=humanAttributes,
       filters="hgnc_symbol", 
       values=human_PPARG_gene, 
       mart=human,
       attributesL=ratAttributes, 
       martL=rat)
colnames(ratMapped_PPARG_Genes) <- c("human_ensemblID", "human_Entrez", "human_Symbol", "rat_ensemblID", "rat_Entrez", "rat_Symbol", "rat_Description", "rat_Chr","rat_Start","rat_End","rat_Strand", "rat_GeneType", "rat_ensemblTranscript")
ratMapped_PPARG_Genes = ratMapped_PPARG_Genes[! (duplicated(ratMapped_PPARG_Genes$rat_ensemblID) & duplicated(ratMapped_PPARG_Genes$human_ensemblID ) ), ]
                  
   
PPAR_genes_extracted = cbind(merged_SstRma[match(ratMapped_PPARG_Genes$rat_Symbol, merged_SstRma[,"Gene Symbol"]),],ratMapped_PPARG_Genes[,c("human_ensemblID","human_Symbol","rat_Symbol","rat_ensemblID","rat_ensemblTranscript")])
PPAR_genes_FAILED_extract = PPAR_genes_extracted[is.na(PPAR_genes_extracted$ID) | PPAR_genes_extracted$"Gene Symbol"=="",]
PPAR_genes_extracted = PPAR_genes_extracted[!is.na(PPAR_genes_extracted$ID) & !PPAR_genes_extracted$"Gene Symbol"=="",]

#write.csv(PPAR_genes_extracted[,c('Gene Symbol',)], file='csvs/check_ortholog_map.csv',row.names=F)

###
write.csv(PPAR_genes_extracted[order(PPAR_genes_extracted$human_Symbol),c('human_Symbol','rat_Symbol','Description','Public Gene IDs','rat_ensemblTranscript')],file='csvs/PPAR_gamma_genes_tested.csv',row.names=FALSE)

###
minPval = matrixStats::rowMins(as.matrix(PPAR_genes_extracted[,c('LA_AvB.P-val','FC_AvB.P-val')]))
minPval_order = order(minPval)
minPval_order = minPval_order[minPval[minPval_order]<0.05]
cols_of_interest=c('ID','rat_Symbol',"human_Symbol",'Public Gene IDs','rat_ensemblTranscript','Description','Chromosome','Strand','Start','Stop','LA_AvB.P-val','LA_BvC.P-val','LA_AvB_log2FC', 'LA_BvC_log2FC','FC_AvB.P-val','FC_BvC.P-val','FC_AvB_log2FC', 'FC_BvC_log2FC')
write.csv(PPAR_genes_extracted[minPval_order,cols_of_interest],file='csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Control_v_Injury_minP_NomSig.csv',row.names=FALSE)

###
minPval = matrixStats::rowMins(as.matrix(PPAR_genes_extracted[,c('LA_BvC.P-val','FC_BvC.P-val')])) 
minPval_order = order(minPval)
minPval_order = minPval_order[minPval[minPval_order]<0.05]
cols_of_interest=c('ID','rat_Symbol',"human_Symbol",'Public Gene IDs','rat_ensemblTranscript','Description','Chromosome','Strand','Start','Stop','LA_AvB.P-val','LA_BvC.P-val','LA_AvB_log2FC', 'LA_BvC_log2FC','FC_AvB.P-val','FC_BvC.P-val','FC_AvB_log2FC', 'FC_BvC_log2FC')
write.csv(PPAR_genes_extracted[minPval_order,cols_of_interest],file='csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Injury_v_InjuryPlusTrt_minP_NomSig.csv',row.names=FALSE)


#### Plot the PPAR-gamma related genes
Cont_v_Injury_genes <- read.csv('csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Control_v_Injury_minP_NomSig.csv',stringsAsFactors = F)
Inj_v_InjPlusTrt_genes <- read.csv('csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Injury_v_InjuryPlusTrt_minP_NomSig.csv',stringsAsFactors = F)
candidate_genes <- unique(c(Cont_v_Injury_genes$rat_Symbol,Inj_v_InjPlusTrt_genes$rat_Symbol))
## Importing normalized data
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

pd$Group = plyr::revalue(pd$Group, c("GroupA"="Control", "GroupB"="Injury", "GroupC"="Injury + Trt"))
dat = cbind(pd, t(rma_signal)  )

## Actual plotting
PPAR_genes_plot = PPAR_genes_extracted[PPAR_genes_extracted$`Gene Symbol` %in% candidate_genes,'ID']

pdf('plots/PPAR_Genes_From_Villapol_2018_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (probe_i in PPAR_genes_plot) {
  
  #Change column name for ggplot2 to work
  ii=match(probe_i, PPAR_genes_extracted$ID)
  
  #  custom_title = paste0( 
  #    "ERC p=",as.character(signif(regionSpecific_mergedStats[ii,'ERC_subset_NoAdj_P.Value'],3)), 
  #    "\n ",fixName ) #custom title
  custom_title=PPAR_genes_extracted[ii,'Gene Symbol']
  
  a = ggplot(dat, aes_string(x = 'Group', y = probe_i, fill='brain_region')) +
    geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
    geom_point(aes(col=`brain_region`),size=5,position = position_jitterdodge(jitter.width=0.3,dodge.width=.85)) + 
    labs(y=paste0(probe_i, " Log2(Signal)"),x='Group', title = custom_title) + 
    theme(legend.position='bottom') +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(size=40,colour='black'))
  print(a)			
}
dev.off()