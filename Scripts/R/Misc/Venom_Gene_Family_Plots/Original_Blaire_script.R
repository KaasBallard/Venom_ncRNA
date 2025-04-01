### Remake the grant figure of SVMP array with ATAC-peaks ###
 
library(BSgenome.Cviridis.custom.CroVir.noSeqNamesMods)
library(tidyverse)
library(rtracklayer)
library(ggcoverage)
library(scales)
library(gggenes)
library(cowplot)
library(viridis)
library(chromVAR)
library(Biostrings)
library(patchwork)
library(scales)
library(ggforce)
library(pheatmap)
 
#### SVMP Array plot (from Blair) ####
# CHANGE PATH UP TO /Dropbox
setwd('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/')
 
pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes,by=c('txid'='X6')) %>%
  mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))
 
 
all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>%
  filter(str_detect(X9,'trnascan',negate = T)) %>%
  mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>%
  mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>%
  filter(tx_id %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>%
  select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>%
  mutate(strand = ifelse(strand == '+','forward','reverse')) %>%
  mutate(direction = ifelse(strand == 'forward',1,-1)) %>%
  mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>%
  left_join(exp) %>%
  mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>%
  mutate(prom_start = ifelse(strand=='forward',start,end))
 
 
 
 
SVMP.info <- all_info %>%
  filter(str_detect(gene,'SVMP|ADAM'))
 
 
 
###
### ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.
###
#SVMP.reg.start <- min(SVMP.info$start)-20000
#SVMP.reg.start <- 13500000
SVMP.reg.start = floor(13901005 / 1000 ) * 1000
#SVMP.reg.start <- 8500000 # SVSP
 
#SVMP.reg.end <- max(SVMP.info$end)+20000
#SVMP.reg.end <- 15000000
SVMP.reg.end = ceiling(14424729 / 1000 ) * 1000
#SVMP.reg.end <- 9000000 # SVSP
 
 
 
 
SVMP.reg.length <- paste(c(round((SVMP.reg.end-SVMP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')
 
# Read in vPERs and super-enhancers
 
svmp.vpers <-
  read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>%
  #read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed', col_names = F) %>%
  select(molecule=1,start=2,end=3,id=4) %>% 
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>%
  mutate(gene = str_split(gene,'\\.')) %>%
  unnest(gene) %>%
  mutate(gene = str_replace(gene,'\\_',' ')) %>%
  left_join(SVMP.info,by='gene') %>%
  mutate(type=' vPERs') %>%
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)
 
#### Get venom gene expression ####
VST_count_mat <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)
 
# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>%
  select(matches('LVG|RVG')) %>%
  filter(grepl('Venom', row.names(.)))
 
 
column_numbers <- 1:13
target_names <- paste0('VG_', column_numbers)
 
for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))
 
  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = TRUE)
    }
  }
}
 
VST_count_mat_filtered <- VST_count_mat %>%
  select(-contains('LVG'),-contains('RVG')) %>%
  select(contains('VG')) %>%
  rownames_to_column(var = 'Gene') %>%
  filter(grepl('SVMP|ADAM', Gene)) %>%
  mutate(AverageExp = rowMeans(select(., starts_with("VG_")), na.rm = TRUE)) %>%
  select(Gene, AverageExp)
 
VST_count_mat_filtered <- VST_count_mat_filtered %>%
  mutate(Gene = gsub('Venom_SVMP', 'SVMP ', Gene)) %>%
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))
 
# ADD TO SVMP INFO
SVMP.info <- SVMP.info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
SVMP.info <- SVMP.info %>% mutate(gene = gsub(' ','',gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ','',gene))
 
SVMP.info$color = "blue"
#### Plotting ####
p1 <- ggplot(SVMP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, color = log10(AverageExp+1))) +
  ggrepel::geom_text_repel(data = all_info %>%
                             filter(str_detect(gene,'SVMP')) %>% # Change for SVSP-SVMP
                             mutate(start = (start + end)/2),
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  #  ggrepel::geom_text_repel(data = all_info %>%
  #                             filter(str_detect(gene,'ADAM')) %>%
  #                             mutate(start = (start + end)/2),
  #                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = -0.5,size=3,color='grey60') +
  geom_segment(aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_diagonal(inherit.aes = F,data=svmp.vpers,aes(x=prom_start,xend=start,y='mi1',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) + # Change for SVSP-SVMP
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  geom_segment(inherit.aes = F,data=svmp.vpers,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svmp.vpers, aes(x=(start+end)/2, y=type),size=2) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
       plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank())