## 05b - Post-hoc PPAR-gamma analysis
# Test an a priori hypothesis about PPAR-gamma genes

setwd('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project')
load('rda/injury_DE_res.rda', verbose=T)
library(dplyr)
library(ggplot2)
library('ggrepel')
theme_set(theme_bw(base_size=24) + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  legend.position="none"				 )) 

#### Genes from Villapol 2018
human_PPARG_gene = c("PPARG","RXRA", "RXRG", "RXRB", "CD36", "JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","FOSL1P1","STAT1", "STAT3", "STAT6", "NFKB1","NFKB2","REL","RELA","RELB", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "MMP1", "MMP3", "MMP9", "IL10", "ARG1", "NOS2","PTGS2","IL6", "TNF","SOD1","SOD2","SOD3", "MAPK1", "MAPK3", "PTGS2")  


### Map human genes found in literature to rat orthologs
library(biomaRt)
ensembl <- useMart("ensembl")
human = useDataset("hsapiens_gene_ensembl", mart=ensembl)
rat = useDataset("rnorvegicus_gene_ensembl", mart=ensembl)


humanAttributes = c("ensembl_gene_id","entrezgene","hgnc_symbol")
ratAttributes = c("ensembl_gene_id","entrezgene","rgd_symbol","description","chromosome_name","start_position","end_position","strand","gene_biotype","ensembl_transcript_id")

ratMapped_PPARG_Genes = getLDS(attributes=humanAttributes,
                               filters="hgnc_symbol", 
                               values=human_PPARG_gene, 
                               mart=human,
                               attributesL=ratAttributes, 
                               martL=rat)
colnames(ratMapped_PPARG_Genes) <- c("human_ensemblID", "human_Entrez", "human_Symbol", "rat_ensemblID", "rat_Entrez", "rat_Symbol", "rat_Description", "rat_Chr","rat_Start","rat_End","rat_Strand", "rat_GeneType", "rat_ensemblTranscript")
#
length(ratMapped_PPARG_Genes$rat_ensemblID)
ratMapped_PPARG_Genes = ratMapped_PPARG_Genes[! (duplicated(ratMapped_PPARG_Genes$rat_ensemblID) & duplicated(ratMapped_PPARG_Genes$human_ensemblID ) ), ] #drop multimapping genes
length(ratMapped_PPARG_Genes$rat_ensemblID)

##
PPAR_genes_extracted = cbind(merged_SstRma_disease[match(as.character(ratMapped_PPARG_Genes$rat_Entrez), as.character(merged_SstRma_disease[,"EntrezID"]) ) ,],ratMapped_PPARG_Genes[,c("human_ensemblID","human_Symbol","rat_Symbol",'rat_Entrez',"rat_ensemblID","rat_ensemblTranscript")])
PPAR_genes_FAILED_extract = PPAR_genes_extracted[is.na(PPAR_genes_extracted$transcriptCluster) | PPAR_genes_extracted$"Symbol"=="" | is.na(PPAR_genes_extracted$rat_Entrez),]
PPAR_genes_extracted = PPAR_genes_extracted[!is.na(PPAR_genes_extracted$transcriptCluster) & !PPAR_genes_extracted$"Symbol"=="" & !is.na(PPAR_genes_extracted$rat_Entrez),]

###
merged_SstRma_disease$PPARg <- ifelse(merged_SstRma_disease$EntrezID %in% PPAR_genes_extracted$rat_Entrez, TRUE,FALSE)

fisher.test(table(merged_SstRma_disease$PPARg,merged_SstRma_disease$Region_MinP<0.0001) )
fisher.test(table(merged_SstRma_disease$PPARg,merged_SstRma_disease$Region_MinP<0.001) )
fisher.test(table(merged_SstRma_disease$PPARg,merged_SstRma_disease$Region_MinP<0.01) )
fisher.test(table(merged_SstRma_disease$PPARg,merged_SstRma_disease$Region_MinP<0.05) )


PPAR_genes_extracted$minP_nomSig = ifelse(PPAR_genes_extracted$Region_MinP<0.05,TRUE,FALSE)
PPAR_genes_extracted$minP_Bonferroni = p.adjust(PPAR_genes_extracted$Region_MinP, method ="bonferroni", n=length(PPAR_genes_extracted$Region_MinP)*2 )
bonf_sig = PPAR_genes_extracted[PPAR_genes_extracted$minP_Bonferroni<0.05,]
bonf_sig[order(bonf_sig$minP_Bonferroni),]

table(PPAR_genes_extracted$minP_nomSig )
table(PPAR_genes_extracted$minP_nomSig, PPAR_genes_extracted$Region_MinP<0.05)
PPAR_genes_extracted[PPAR_genes_extracted$minP_nomSig,]

cols_of_interest=c('transcriptCluster','Symbol',"human_Symbol",'Public Gene IDs','rat_ensemblTranscript','Description','Chromosome','Strand','Start','Stop','LA_AvB.P-val','LA_BvC.P-val','LA_AvB_log2FC', 'LA_BvC_log2FC','FC_AvB.P-val','FC_BvC.P-val','FC_AvB_log2FC', 'FC_BvC_log2FC')
write.csv(PPAR_genes_extracted[minPval_order,cols_of_interest],file='csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Control_v_Injury_minP_NomSig.csv',row.names=FALSE)

#### Plot the PPAR-gamma related genes
Cont_v_Injury_genes <- read.csv('csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Control_v_Injury_minP_NomSig.csv',stringsAsFactors = F)
Inj_v_InjPlusTrt_genes <- read.csv('csvs/PPAR_gamma_Villapol2018_rat_gene_statistics_sortedBy_Injury_v_InjuryPlusTrt_minP_NomSig.csv',stringsAsFactors = F)
candidate_genes <- unique(c(Cont_v_Injury_genes$rat_Symbol,Inj_v_InjPlusTrt_genes$rat_Symbol))
## Importing normalized data
RMA_Res = list.files('C:/Users/steph/Documents/R/UMB/PTSD_Rat_Project/TAC/SST-RMA_BothRegions',pattern = '.txt',full.names = T)

res = lapply(RMA_Res, function(x) as.data.frame(data.table::fread(x)))
names(res) = gsub("\\..*", "", basename(RMA_Res))


## Actual plotting
PPAR_genes_plot = PPAR_genes_extracted[PPAR_genes_extracted$`Gene Symbol` %in% candidate_genes,'ID']

pdf('plots/PPAR_Genes_From_Villapol_2018_boxplots.pdf',height=10,width=12,useDingbats=FALSE)
for (probe_i in PPAR_genes_plot) {
  
  #Change column name for ggplot2 to work
  ii=match(probe_i, PPAR_genes_extracted$ID)
  
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