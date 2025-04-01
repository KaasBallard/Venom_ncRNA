### Remake the grant figure of PLA2 array with ATAC-peaks ###

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

#### PLA2 Array plot (from Blair) ####
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




PLA2.info <- all_info %>% 
  filter(str_detect(gene,'PLA2'))



###
### ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions. 
###
#SVMP.reg.start <- min(SVMP.info$start)-20000
#SVMP.reg.start <- 13500000
#PLA2.reg.start = floor(3019890 / 1000 ) * 1000
PLA2.reg.start = 3023000
#SVMP.reg.start <- 8500000 # SVSP

#SVMP.reg.end <- max(SVMP.info$end)+20000
#SVMP.reg.end <- 15000000
PLA2.reg.end = ceiling(3044567 / 1000 ) * 1000 
#SVMP.reg.end <- 9000000 # SVSP




PLA2.reg.length <- paste(c(round((PLA2.reg.end-PLA2.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in vPERs and super-enhancers

PLA2.vpers <- 
  read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/PLA2_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  #read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  #read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed', col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(PLA2.info,by='gene') %>% 
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
  filter(grepl('PLA2', Gene)) %>% 
  mutate(AverageExp = rowMeans(select(., starts_with("VG_")), na.rm = TRUE)) %>% 
  select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  mutate(Gene = gsub('Venom_PLA2', 'PLA2 ', Gene)) %>% 
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))

# ADD TO SVMP INFO
PLA2.info <- PLA2.info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
PLA2.info <- PLA2.info %>% mutate(gene = gsub(' ','',gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ','',gene))

CTCF_motifs <- read.table('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/analysis/10_loops_chromatin/contact_loops/CTCF_BoundMotifs_02.03.21.bed', header = F, col.names = c('chr', 'start','end','name'))
CTCF_motifs <- CTCF_motifs %>% 
  filter(chr == 'scaffold-mi7') %>% 
  mutate(site = (start + end)/2) %>% 
  filter(start < PLA2.reg.end & start > PLA2.reg.start) %>% 
  mutate(type = '  CTCF') # weird fix to change the order of plotting????



#### Plotting ####
p1 <- ggplot(PLA2.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(AverageExp+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'PLA2')) %>% # Change for SVSP-SVMP-PLA2
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 1,size=3) +
  #  ggrepel::geom_text_repel(data = all_info %>% 
  #                             filter(str_detect(gene,'ADAM')) %>% 
  #                             mutate(start = (start + end)/2), 
  #                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = -0.5,size=3,color='grey60') +
  geom_segment(aes(x=PLA2.reg.start,xend=PLA2.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_diagonal(inherit.aes = F,data=PLA2.vpers,aes(x=prom_start,xend=start,y='mi7',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) + # Change for SVSP-SVMP
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  geom_segment(inherit.aes = F,data=PLA2.vpers,aes(x=PLA2.reg.start,xend=PLA2.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=PLA2.vpers, aes(x=(start+end)/2, y=type),size=2) +
  #geom_segment(inherit.aes = F,data=CTCF_motifs,aes(x=PLA2.reg.start,xend=PLA2.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  #geom_point(inherit.aes = F, data=CTCF_motifs, aes(x=site, y=type),size=2, color = 'darkorchid3') +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(PLA2.reg.start,PLA2.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank())


#### PLA2 Array ATAC-seq with peaks ####
setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/6_GeneVignettes/PLA2_wholeArray')
peakfile <- paste('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/merged_peaks/allSample_mergedPeaks.bed')
peaks <- getPeaks(peakfile, sort_peaks = FALSE)
peaks <- resize(peaks, width = 500, fix = "center")
peaks <- as.data.frame(peaks) %>% filter(seqnames == 'scaffold-mi7')

# Load data and plot
track.folder = '/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/bams'

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder), value = T)))
sample.meta <- sample.meta %>% 
  mutate(Type = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Group = str_extract(SampleName, 'viridis|concolor|lutosus|cerberus')) %>% 
  filter(Group != 'lutosus') %>% # remove lutosus for now because its ATAC-seq quality is bad
  #filter(Type != 'CV1096') %>% # remove CV1096 for now because we don't have a genome
  arrange(match(Group, c('viridis', 'concolor', 'cerberus')))

# load regions of interest from bam files
chrom = "scaffold-mi7"
start_pos = 3023000
end_pos = ceiling(3044567 / 1000 ) * 1000
track.df = LoadTrackFile(track.folder = track.folder, format = "bam", # change to bam
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = T, # change to T for BAM
                         single.nuc.region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change region to single.nuc.region for BAM, change to region for bw
)

score_variance <- track.df %>%
  group_by(seqnames, start, end) %>%
  summarise(score = var(score))
score_variance$Type <- 'Variance'
score_variance$Group <- 'Variance'


# Add BAM info # For genome bams only
track.folder_genomes = '/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/genomic_bams/' # change from bams (ATAC) to genomic bams
sample.meta_genomes = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder_genomes), value = T)))
sample.meta_genomes <- sample.meta_genomes %>% 
  mutate(Type = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Group = str_extract(SampleName, 'viridis|concolor|lutosus|cerberus')) %>% 
  filter(Group != 'lutosus') %>% # remove lutosus for now because its ATAC-seq quality is bad
  #filter(Type != 'CV1096') %>% # remove CV1096 for now because we don't have a genome
  arrange(match(Group, c('viridis', 'concolor', 'cerberus')))

average_coverage <- data.frame(
  stringsAsFactors = FALSE,
  Sample_Name = c("CV0857","CV0985","CV0987",
                  "CV1096","CV1084","CV1081","CV1089","CV1090","CV1087",
                  "CV1086","CV1082","CV1095"),
  Average.coverage = c(65.74,52.18,55.56,57.00,47.7,
                       58.01,7.73,54.36,47.77,42.1,54.01,59.95)
)


sample.meta_genomes <- sample.meta_genomes %>% left_join(average_coverage, by = c('Type'='Sample_Name'))
target_coverage <- mean(sample.meta_genomes$Average.coverage)
scaling_factors <- target_coverage / sample.meta_genomes$Average.coverage
sample.meta_genomes$Coverage.Scaling <- scaling_factors

track.df_genome = LoadTrackFile(track.folder = track.folder_genomes, format = "bam", # change to bam
                                bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                                meta.info = sample.meta_genomes,
                                single.nuc = T, # change to T for BAM
                                single.nuc.region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change region to single.nuc.region for BAM, change to region for bw
)
# Modify genome coverage scores so that each genome effectively has the same depth (otherwise low depth genome will look like they have no coverage, when really they do, it's just a library size difference)
track.df_genome$scaled_score <- track.df_genome$score * track.df_genome$Coverage.Scaling

# Add columns for data type
track.df$Data_Type <- "ATAC-seq"
track.df_genome$Data_Type <- "genomic"

combined_df <- bind_rows(track.df, track.df_genome)
combined_df <- combined_df %>% 
  mutate(final_score = case_when(!is.na(score) ~ score,
                                 !is.na(scaled_score) ~ scaled_score,
                                 TRUE~NA)) # move the genomic and atac scores to a column to be plotted


p3 <- 
  ggplot() +
  theme_minimal() +
  #geom_col(data = track.df, aes(x = start, y = score, fill = Group)) + 
  geom_col(data = combined_df %>% filter(Group == 'concolor'), aes(x = start, y = final_score, fill = Data_Type)) + # keep only concolor for now
  #facet_wrap(~factor(Type, levels = c(sample.meta$Type, 'Variance')), ncol = 1, strip.position = "right", scales = "free_y") + # see all samples
  #facet_wrap(~factor(Data_Type), ncol = 1, scales = 'free_y') + # facet by genomic vs atac
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        #strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, title = NULL)) +
  # scale_fill_manual(values = c("viridis" = "black", # #2E8B58
  #                              "concolor" = "black", # #D9A528
  #                              "cerberus" = "black",
  #                              "Variance" = "blue"),
  #                   labels = c(expression(italic("C. viridis")),
  #                              expression(italic("C. concolor")),
  #                              expression(italic("C. cerberus")),
  #                              breaks = c("viridis", "concolor", "cerberus"))) +
  scale_fill_manual(values = c("ATAC-seq" = "black",
                               "genomic" = "red")) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "read density") +
  ggtitle("Accessibility at PLA2 array")
x <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files/Cvv_GTF_to_converted_names_05.22.23.txt', header = T)
# Plot variance
p4 <- ggplot() +
  theme_minimal() +
  geom_col(data = score_variance, aes(x = start, y = score), fill = 'blue') + 
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, title = NULL))

plot_grid(p1, p2, p3, p4, nrow = 4, axis = 'lr', align = 'hv', rel_heights = c(0.4,1,0.6,0.2))

#### Make Individual Gene Vignettes ####
setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/6_GeneVignettes/')

peakfile <- paste('../1_chromVAR/merged_peaks/allSample_mergedPeaks.bed', sep = '/')
peaks <- getPeaks(peakfile, sort_peaks = FALSE)
peaks <- resize(peaks, width = 500, fix = "center")

# Preform overlap with venom enhancers and promoters
enh_pr <- read.table('../../z_data_files/venom_enhancers_promoters-1000+100_09.12.23.txt', header = T) # change this
enh_pr_gr <- with(enh_pr, GRanges(seqnames, IRanges(start, end))) # the original file is -500bp + 10bp around the venom exon 1, make it -1000bp +100bp
overlap <- findOverlaps(enh_pr_gr, peaks)
enh_pr_peaks <- cbind(as.data.frame(enh_pr[queryHits(overlap),]), peaks[subjectHits(overlap),])
colnames(enh_pr_peaks) <- c("seqnames", "feature_start", "feature_end", "name", "type", "seqnames2", "peak_start", "peak_end", "width", "strand")
enh_pr_peaks <- enh_pr_peaks %>% dplyr::select(-seqnames2, strand)
rm(enh_pr_gr, overlap) # remove temp vars

# Load data and plot
track.folder = '/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/bams/' # change from bams (ATAC) to genomic bams

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder), value = T)))

# all samples
sample.meta <- sample.meta %>% 
  mutate(Type = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Group = str_extract(SampleName, 'viridis|concolor|lutosus|cerberus')) %>% 
  filter(Group != 'lutosus') %>% # remove lutosus for now because its ATAC-seq quality is bad
  #filter(Type != 'CV1096') %>% # remove CV1096 for now because we don't have a genome
  arrange(match(Group, c('viridis', 'concolor', 'cerberus')))


# load regions of PLA2 interest from bam files
feature_of_interest = 'PER39_PLA2A1'

chrom = "scaffold-mi7"
start_pos = floor(enh_pr[enh_pr$name==feature_of_interest,'start'] / 500 ) * 500 
end_pos = ceiling(enh_pr[enh_pr$name==feature_of_interest,'end'] / 500 ) * 500 

#start_pos = 
end_pos = end_pos + 500

track.df = LoadTrackFile(track.folder = track.folder, format = "bam", # change to bam
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = T, # change to T for BAM
                         single.nuc.region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change region to single.nuc.region for BAM,
)

#TFBS_table_e <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/2_TFBS_Presence_Absence/z_outputs/All_SNPs/PER_ciiider_table_allSNPs_fp_SNP_PA_Sept2023.txt.gz', header = T) # change from promoter to PER
#TFBS_table_p <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/2_TFBS_Presence_Absence/z_outputs/All_SNPs/Promoter_ciiider_table_allSNPs_fp_SNP_PA_Aug2023.txt.gz', header = T) # change from promoter to PER
TFBS_table_f <- TFBS_table_p %>% 
  filter(gene == feature_of_interest) %>% # change from Promoter to PER18
  select(Sample.Name, Transcription.Factor.Name, genome.TFBS.start, genome.TFBS.end) %>% 
  mutate(Sample.Name = gsub("_.*", "", Sample.Name)) %>% 
  dplyr::rename(Type = Sample.Name) %>% 
  #mutate(TF.Name2 = ifelse(Transcription.Factor.Name %in% c('NFATC1', 'NFIX', 'ZBTB26'), Transcription.Factor.Name, 'Other')) %>% # I think NFATC1/NFIX candidates in promoter
  #mutate(TF.Name2 = ifelse(Transcription.Factor.Name %in% c("Arnt", "CLOCK", "Creb3l2", "EHF", "ELF5", "ELK4", "FIGLA", "GABPA", "HES6", "Irx2", "KLF11", "KLF16", "MEIS1", "NFATC1", "NFIA", "NR1H4", "RORC"), Transcription.Factor.Name, 'Other')) %>% # I think these candidates in enhancer
  #arrange(factor(TF.Name2, levels = c("Other", unique(Transcription.Factor.Name)))) %>% 
  distinct()
custom_colors <- c("Other" = "grey60", 
                   "NFATC1" = "#E66100", 
                   "NFIX" = "#E1BE6A",
                   "ZBTB26" = "#4992C2") # promoter
custom_colors <- c("Arnt" = "#A6CEE3", "CLOCK" = "#4992C2", "Creb3l2" = "#569EA4",
                   "EHF" = "#AADB84", "ELF5" = "#52AF43", "ELK4" = "#8A9D5B",
                   "FIGLA" = "#F88A89", "GABPA" = "#E73233", "HES6" = "#F06C45",
                   "Irx2" = "#FDB35A", "KLF11" = "#FE870D", "KLF16" = "#E19B78",
                   "MEIS1" = "#B294C7", "NFATC1" = "#70449D", "NFIA" = "#C7B699",
                   "NR1H4" = "#E6CB75", "RORC" = "#B15928", "Other" = "grey60") # PER17/PER18


peak_shift <- 0
rect <- peaks %>% 
  as.data.frame() %>% 
  filter(seqnames == 'scaffold-mi7') %>% 
  select(start, end) %>% 
  mutate(peak_start = start + peak_shift, peak_end = end + peak_shift) %>% 
  mutate(Type = 'CV1090')

example_genes <- enh_pr %>%
  filter(name == feature_of_interest) %>%
  select(start, end, name) %>%
  mutate(Type = 'CV1090')

p3 <- 
  ggplot() +
  theme_minimal() +
  geom_col(data = track.df, aes(x = start, y = score, fill = Group)) +
  annotate(geom = "rect",
           xmin = rect$peak_start,
           xmax = rect$peak_end,
           ymin = -Inf,
           ymax = +Inf,
           alpha = 0.2) +
  geom_point(data = TFBS_table_f,
             aes(x = genome.TFBS.start, y = 100, fill = Transcription.Factor.Name, color = Transcription.Factor.Name),
             size = 2,
             position = position_jitter(width=0, height=10),
             inherit.aes = F) +
  #scale_color_manual(values = custom_colors) +
  geom_gene_arrow(data = example_genes, inherit.aes = F, aes(xmin = start, xmax = end, y = 40, fill="blue"), arrowhead_height = grid::unit(2, "mm"), arrow_body_height = grid::unit(1, "mm")) +
  facet_wrap(~factor(Type, levels = sample.meta$Type), ncol = 1, strip.position = "right") +
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6)) +
  scale_y_continuous(breaks = c(0,200)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        #strip.text.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, title = NULL)) +
  scale_fill_manual(values = c("viridis" = "#2E8B58",
                               "concolor" = "#D9A528",
                               "cerberus" = "black"),
                    labels = c(expression(italic("C. viridis")),
                               expression(italic("C. concolor")),
                               expression(italic("C. cerberus")),
                               breaks = c("viridis", "concolor", "cerberus"))) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq read density") +
  #ggtitle(paste0(rect$TF, " associated variability")) 
  ggtitle(paste0("Accessibility landscape\nat ",feature_of_interest))


#### Gene Expression ####
VST_count_mat <- read.table('../../z_data_files/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T) # use either the VST normalized mat or the normalized count mat

# Log Normalize counts if using normalized counts
# VST_count_mat =  log(VST_count_mat + 1)

# Filter for Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
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
  select(-contains('LVG'),-contains('RVG'))

samples <- read.table('../../z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
samples <- samples %>% 
  mutate(rna_12_id = str_extract(Sample_ID, str_extract(Sample_ID, "\\d+")))

other_meta <- readxl::read_xlsx('../../z_data_files/fixed_12snake_metadata_05.21.23.xlsx', sheet = 2, col_names = T) %>% 
  select(Corrected_name, 7) %>% 
  dplyr::slice(1:13) %>% 
  janitor::clean_names() %>% 
  dplyr::slice(-10)


# all samples
samples <- samples %>% 
  left_join(other_meta, by = 'rna_12_id') %>% 
  filter(grepl('RVG', Sample_ID)) %>% 
  select(Lineage, Population, corrected_name, rna_12_id) %>% 
  dplyr::rename(SampleID = corrected_name) %>% 
  relocate(SampleID, .before=Lineage) %>% 
  filter(!grepl('CV1084', SampleID)) %>% # remove bad viridis sample
  filter(Lineage != 'lutosus') %>% # remove lutosus with no ATAC
  mutate(rna_12_id = paste0('VG_', rna_12_id)) %>% 
  arrange(match(SampleID, sample.meta$Type)) # match the order of the ATAC-seq tracks

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  select(any_of(samples$rna_12_id))
colnames(VST_count_mat_filtered) <- samples$SampleID[match(colnames(VST_count_mat_filtered), samples$rna_12_id)]
VST_count_mat_filtered %>% 
  select(one_of(samples$SampleID), everything())

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  filter(rownames(.) == GEX_of_interest) %>% 
  t() %>% 
  as.data.frame()
VST_count_mat_filtered$SampleID <- rownames(VST_count_mat_filtered)
p4 <- VST_count_mat_filtered %>% 
  left_join(samples %>% select(SampleID, Population, Lineage), by = 'SampleID') %>% 
  ggplot(., aes(x = SampleID, y = Venom_PLA2A1, fill = Lineage)) +
  theme_bw() +
  geom_col() +
  labs(y = "normalized\ncounts", x = NULL) +
  coord_flip() +
  scale_x_discrete(limits = rev(samples$SampleID),
                   labels = rev(samples$SampleID)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = c("viridis" = "#2E8B58",
                               "concolor" = "#D9A528",
                               "cerberus" = "black")) +
  guides(fill = 'none')


plot_grid(p3, p4, nrow = 1, align = 'hv', axis = 'tb', rel_widths = c(3,1))

