#### Set up and Read Data ####

# Load in packages
library(cowplot)
library(tidyverse)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)
library(ggplot2)

# Set working directory
setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
#setwd('C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA')

# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.tsv'
# miRNA_mRNA_data <- 'Data/Merged/miRNA_mRNA_Combined_Data_IMPORTANT.tsv'

# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)
# miRNA_mRNA_data <- read.table(file = miRNA_mRNA_data, header = T)

# Create shorter df name
mi_df <- miRNA_mRNA_protein_df %>% 
  rename('Genes' = 'Converted.Gene.IDs')


#### PLA2 Array plot ####

# Load packages
library(tidyverse)
library(rtracklayer)
library(ggcoverage)
library(scales)
library(gggenes)
library(viridis)
library(ggforce)


# CHANGE PATH UP TO /Dropbox
setwd('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/')

# Read priority venom genes data
pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)

# Read gene expression data and filter for priority venom genes
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes,by=c('txid'='X6')) %>%
  mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))

# Read all gene information and filter for priority venom genes
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

# Filter for PLA2 genes
PLA2_info <- all_info %>%
  filter(str_detect(gene,'PLA2'))


# ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.

# Define region start and end points for PLA2 plot
PLA2_reg_start = floor(3019890 / 1000 ) * 1000
PLA2_reg_end = ceiling(3043778 / 1000 ) * 1000
PLA2_reg_length <- paste(c(round((PLA2_reg_end-PLA2_reg_start)/1000,digits = 2),'kb'),collapse = ' ')


# Get venom gene expression 
VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Data/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)

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
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
    }
  }
}

VST_count_mat_filtered <- VST_count_mat %>%
  select(-contains('LVG'),-contains('RVG')) %>%
  select(contains('VG')) %>%
  #select('RVG_5', 'RVG_6', 'RVG_7', 'RVG_12') %>%
  rownames_to_column(var = 'Gene') %>%
  filter(grepl('PLA2', Gene)) %>%
  mutate(AverageExp = rowMeans(select(., starts_with("VG_")), na.rm = T)) %>%
  select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>%
  mutate(Gene = gsub('Venom_PLA2', 'PLA2 ', Gene)) %>%
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))

# ADD TO PLA2 INFO
PLA2_info <- PLA2_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
PLA2_info <- PLA2_info %>% mutate(gene = gsub(' ','',gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ','',gene))

# Add venom gene expression data to PLA2_info
PLA2_info <- PLA2_info %>%
  mutate(
    # Swap start and end if direction is -1
    start2 = case_when(direction == -1 ~ end,
                       T ~ start),
    end2 = case_when(direction == -1 ~ start,
                     T ~ end)
  )

# Set color for PLA2_info
PLA2_info$color = "blue"

# Plotting 
# Create plot using ggplot
PLA2_figure <- ggplot(PLA2_info, aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp+1))) +
  # Add text labels for genes
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% 
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2,size=2.5) +
  # Add gray segment indicating the region limits                         
  geom_segment(aes(x=PLA2_reg_start,xend=PLA2_reg_end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  # Add gene arrows
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  # Set labels and scales
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(PLA2_reg_start,PLA2_reg_end),expand=c(0,0)) +
  # Set theme
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank())

# Save plot to a PDF file
ggsave("/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Scripts/R/SVMP_plot/PLA2_array_plot.pdf", PLA2_figure, width = 10, height = 2, units = "in")



#### PLA2 Array plot with miRNA arrows ####

# Load packages
library(tidyverse)
library(rtracklayer)
library(ggcoverage)
library(scales)
library(gggenes)
library(viridis)
library(ggforce)

# Set up the first part of the figure for the gene arrows
# CHANGE PATH UP TO /Dropbox
setwd('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/')

# Read priority venom genes data
pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)

# Read gene expression data and filter for priority venom genes
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes,by=c('txid'='X6')) %>%
  mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))

# Read all gene information and filter for priority venom genes
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

# Filter for PLA2 genes
PLA2_info <- all_info %>%
  filter(str_detect(gene,'PLA2'))


# ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.

# Define region start and end points for PLA2 plot
PLA2_reg_start = floor(3019890 / 1000 ) * 1000
PLA2_reg_end = ceiling(3043778 / 1000 ) * 1000
PLA2_reg_length <- paste(c(round((PLA2_reg_end-PLA2_reg_start)/1000,digits = 2),'kb'),collapse = ' ')


# Get venom gene expression 
VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)

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
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
    }
  }
}

VST_count_mat_filtered <- VST_count_mat %>%
  select(-contains('LVG'),-contains('RVG')) %>%
  select(contains('VG')) %>%
  #select('RVG_5', 'RVG_6', 'RVG_7', 'RVG_12') %>%
  rownames_to_column(var = 'Gene') %>%
  filter(grepl('PLA2', Gene)) %>%
  mutate(AverageExp = rowMeans(select(., starts_with("VG_")), na.rm = T)) %>%
  select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>%
  mutate(Gene = gsub('Venom_PLA2', 'PLA2 ', Gene)) %>%
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))

# ADD TO PLA2 INFO
PLA2_info <- PLA2_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
PLA2_info <- PLA2_info %>% mutate(gene = gsub(' ','',gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ','',gene))

# Add venom gene expression data to PLA2_info
PLA2_info <- PLA2_info %>%
  mutate(
    # Swap start and end if direction is -1
    start2 = case_when(direction == -1 ~ end,
                       T ~ start),
    end2 = case_when(direction == -1 ~ start,
                     T ~ end)
  )

# Set color for PLA2_info
PLA2_info$color = "blue"


# Set stuff up for the miRNA portion of the figure
# Create new data frame with only miRNA count data.
PLA2_miRNAs_df <- mi_df %>% 
  filter(!Origin %in% c('exon', 'five_prime_utr')) %>%  # Filter out anything that didn't come from the 3'UTR
  filter(!miRNA.Counts.RVG.5S.CV1087.viridis.North.F == 0) %>% 
  select(miRNA.Cluster, Genes, miRNA.Counts.RVG.5S.CV1087.viridis.North.F) %>% 
  rename_with(~ "counts", 'miRNA.Counts.RVG.5S.CV1087.viridis.North.F') %>% 
  #rename(counts = 'miRNA.Counts.RVG.5S.CV1087.viridis.North.F') # I have no idea why this doesn't work but it doesn't
  mutate(Genes = str_replace(Genes, 'Venom_', '')) %>% # Remove the venom wording so that it matches what is in PLA2_info
  mutate(Genes = str_replace(Genes, 'PLA2G2E.1', 'PLA2gIIE')) %>% 
  distinct()


# Merge PLA2_info to SVSP-miRNAs_df
PLA2_miRNA_info <- left_join(PLA2_info, PLA2_miRNAs_df, by = c('gene' = 'Genes'))

# Create new data frame so I can do this in illustrator :(
filt_PLA2_info <- PLA2_miRNA_info %>% 
  select(gene, miRNA.Cluster)

# Plotting 
# Create plot using ggplot
PLA2_figure <- ggplot(PLA2_info, aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp+1))) +
  # Add text labels for genes
  ggrepel::geom_text_repel(data = all_info %>%
                             filter(str_detect(gene,'PLA2')) %>%
                             mutate(start = (start + end)/2),
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2, size=2.5) +
  ggrepel::geom_text_repel(data = PLA2_miRNA_info,
                           aes(x = (start + end)/2, label = miRNA.Cluster), nudge_y = -0.5, size = 2.5) +
  # Add gray segment indicating the region limits                         
  geom_segment(aes(x=PLA2_reg_start, xend=PLA2_reg_end, y=str_remove(molecule, 'scaffold\\-'), yend=str_remove(molecule, 'scaffold\\-')), lwd=1, color='grey70') +
  # Add gene arrows
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  # Set labels and scales
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma, limits=c(PLA2_reg_start, PLA2_reg_end), expand=c(0,0)) +
  # Set theme
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black', face='bold', size = 14),
        axis.title.x=element_blank())


# # Plotting 
# # Create plot using ggplot
# PLA2_figure <- ggplot(PLA2_info, aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp+1))) +
#   # Add text labels for genes
#   ggrepel::geom_text_repel(data = all_info %>%
#                              filter(str_detect(gene,'PLA2')) %>%
#                              mutate(start = (start + end)/2),
#                            aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2, size=2.5) +
#   # Add gray segment indicating the region limits
#   geom_segment(aes(x=PLA2_reg_start, xend=PLA2_reg_end, y=str_remove(molecule, 'scaffold\\-'), yend=str_remove(molecule, 'scaffold\\-')), lwd=1, color='grey70') +
#   # Add gene arrows
#   geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
#   # Set labels and scales
#   ylab('') +
#   xlab('') +
#   scale_fill_viridis_c(option = 'B') +
#   scale_x_continuous(labels = comma, limits=c(PLA2_reg_start, PLA2_reg_end), expand=c(0,0)) +
#   # Set theme and remove the scaling in the y-axis
#   theme_classic(base_size = 14) +
#   theme(axis.line.y = element_blank(),
#         plot.title.position = 'plot',
#         plot.title = element_text(color='black', face='bold', size = 14),
#         axis.title.x=element_blank()) 
# # the y-axis is alphabetically ordered


# Save plot to a PDF file
ggsave("/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Figures/Gene_Arrays/PLA2_array_plot_cartoon.pdf", PLA2_figure, width = 10, height = 2, units = "in")
