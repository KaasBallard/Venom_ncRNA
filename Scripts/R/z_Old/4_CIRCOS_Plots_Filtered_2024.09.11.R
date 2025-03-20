# Last Edited: 2024/9/11

#### Set up and Read Data ####

# Load in packages
library(compositions)
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
# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")
# setwd('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')


# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.08.31.tsv'

# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)
# # I am going to exclude the following genes because they are not well annotated
# excluded_genes = c(
#   'maker-scaffold-mi1-augustus-gene-59.13_crovir-transcript-12940',
#   'maker-scaffold-mi1-augustus-gene-59.20_crovir-transcript-12947',
#   'maker-scaffold-mi2-augustus-gene-22.17_crovir-transcript-739',
#   'maker-scaffold-un11-augustus-gene-5.19',
#   'XP_016876419',
#   'XP_011528471'
# )
# # I removed the above method because it isn't capturing every instance of the problem.


# Let's keep the analysis limited to the 3UTR
# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  dplyr::rename(
    'miRNA.Cluster.Original' = 'miRNA.Cluster',
    'Genes' = 'Converted.Gene.IDs',
    'miRNA.Cluster' = 'Putative.miRNA.Name'
  ) %>%
  dplyr::select(miRNA.Cluster, everything()) %>% # Move the new miRNA.Clusters to the front
  filter(!str_detect(Genes, 'maker-scaffold|augustus|XP_')) %>% # This should get rid of all of the weirdly annotated genes
  # filter(!(Genes %in% excluded_genes)) %>% 
  # filter(str_detect(Origin, 'three_prime_utr|five_prime_utr')) %>% # Filter out CDS targeting miRNAs
  dplyr::select(-'gtf.gene', -'crovir.transcript', -'Protein.Probability', -'Top.Peptide.Probability', -'Blast.Alignment.Length')  # Remove columns to save memory
rm(miRNA_mRNA_protein_df)

# Add more filters the main data frame separately so things don't get to complicated
mi_df <- mi_df %>% 
  filter(Total.Energy < -19) %>% # I am including this filter to remove any low strength binding examples
  filter(Total.Score > 150) %>% # I am including this filter to remove any low confidence binding examples
  mutate(
    Line.Type = case_when(
      Origin == 'three_prime_utr' ~ 'solid',
      Origin == 'five_prime_utr' ~ 'longdash',
      Origin == 'CDS' ~ 'twodash',
    )
  ) 

# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
ADAM_color <- '#9A70B5'
SVSP_color <- '#F0B830' 
PLA2_color <- '#7570B3'
VEGF_color <- '#74ADD1'
ohanin_color <- '#3A489C'
myotoxin_color <- '#B2182B'
vQC_color <- '#80BC50'
CRISP_color <- '#E7298A'
CTL_color <- '#F67E17'
EXO_color <- '#49FFFF'
LAAO_color <- '#B35806'
BPP_color <- '#1B9E77'
other_color <- '#666666'
miRNA_color <- 'grey'
three_prime_color <- 'black'
# five_prime_color <- '#0072b2'
five_prime_color <- '#1B9E77'
# cds_color <- '#d55e00'
cds_color <- '#4A70B5'




#### Create General Dataframe for CIRCOS ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/CIRCOS_Functions.R')

# Create smaller data frame for circos plot circos.link function to draw lines
venom_circos_df <- mi_df %>% 
  dplyr::select(
    miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, Positions, Origin, Line.Type,
    contains('miRNA.RPM.'),
    contains('RNA.'),
    contains('Intensity.')
  ) %>% 
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  filter(!miRNA.Target.Chrom %in% c('PE-reconstructed-10x-myo')) %>% 
  filter(!str_detect(miRNA.Sequence.Chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  filter(!str_detect(miRNA.Target.Chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  mutate(Positions = strsplit(as.character(Positions), ' ')) %>% # Separate the Positions column into multiple rows
  unnest(Positions) %>% 
  mutate(Positions = as.numeric(Positions)) %>%  # Convert Positions to numeric
  filter(!is.na(Positions)) %>% 
  # Create the new column miRNA.Target.Position
  mutate(miRNA.Target.Position = miRNA.Target.Start + Positions) %>% 
  mutate(Color = case_when(
    str_starts(Genes, 'Venom_PLA2') ~ PLA2_color,
    str_starts(Genes, 'Venom_SVSP') ~ SVSP_color,
    str_starts(Genes, 'Venom_SVMP') ~ SVMP_color,
    str_starts(Genes, 'Venom_ADAM') ~ ADAM_color,
    str_starts(Genes, 'Venom_vQC') ~ vQC_color,
    str_starts(Genes, 'Venom_ohanin') ~ ohanin_color,
    str_starts(Genes, 'Venom_LAAO') ~ LAAO_color,
    str_starts(Genes, 'Venom_CTL') ~ CTL_color,
    str_starts(Genes, 'Venom_CRISP') ~ CRISP_color,
    str_starts(Genes, 'Venom_BPP') ~ BPP_color,
    str_starts(Genes, 'Venom_VEGF') ~ VEGF_color,
    str_starts(Genes, 'Venom_myotoxin') ~ myotoxin_color,
    str_starts(Genes, 'Venom_EXO') ~ EXO_color,
    T ~ other_color
  ))

# Set track height
trackheight <- 0.2

#### CIRCOS Plot for RVG_5S colored by venom gene family ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/CIRCOS_Functions.R')


# Create smaller data frame for circos plot circos.link function to draw lines
# This data frame also only has RVG_5S information in it and is color coded for venom gene families
RVG5_venom_circos_df <- venom_circos_df %>% 
  dplyr::select(
    -contains('.RVG.6S.'),
    -contains('.RVG.7S.'),
    -contains('.LVG.2.'),
    -contains('.LVG.4.'),
    -contains('.LVG.9.')
  ) %>% 
  filter(!miRNA.RPM.RVG.5S.CV1087.viridis.North.F == 0) %>% # filter out any miRNA with zero expression, I don't think there are any, but just to be safe
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  distinct()

# Initialize max count values
max_protein_value <- range(log(RVG5_venom_circos_df$Intensity.RVG.5S.CV1087.viridis.North.F + 1), na.rm = T)
max_miRNA_count <- range(log(RVG5_venom_circos_df$miRNA.RPM.RVG.5S.CV1087.viridis.North.F + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG5_venom_circos_df$RNA.VST.RVG.5S.CV1087.viridis.North.F + 1), na.rm = T)

# Define breaks for color scales
protein_breaks <- quantile(log(RVG5_venom_circos_df$Intensity.RVG.5S.CV1087.viridis.North.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
miRNA_breaks <- quantile(log(RVG5_venom_circos_df$miRNA.RPM.RVG.5S.CV1087.viridis.North.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
mRNA_breaks <- quantile(log(RVG5_venom_circos_df$RNA.VST.RVG.5S.CV1087.viridis.North.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)

# Create color ramps using magma color map
protein_color_ramp <- colorRamp2(protein_breaks, magma(3))
miRNA_color_ramp <- colorRamp2(miRNA_breaks, magma(3))
mRNA_color_ramp <- colorRamp2(mRNA_breaks, magma(3))



# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/Filtered/Venom_Genes_Colored_CIRCOS_CV1087_viridis_2024.09.02.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par('track.height' = trackheight,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = c(rep(1, 17), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# # Draw circos tracks for protein
# draw_circos_tracks_points2(
#   RVG5_venom_circos_df,
#   max_value = max_protein_value,
#   track_height = 0.2,
#   chromosome_col = 'miRNA.Target.Chrom',
#   loci_location_col =  'miRNA.Target.Start',
#   value_col = 'Intensity.RVG.5S.CV1087.viridis.North.F',
#   color = 'purple'
#   )

# Draw circos tracks for mRNA
draw_circos_tracks_points2(
  RVG5_venom_circos_df,
  max_value = max_mRNA_expression,
  track_height = trackheight,
  chromosome_col = 'miRNA.Target.Chrom',
  loci_location_col = 'miRNA.Target.Start',
  value_col = 'RNA.VST.RVG.5S.CV1087.viridis.North.F',
  color_ramp = mRNA_color_ramp,
  chrom_labels = T
)

# Draw circos tracks for miRNA
draw_circos_tracks_points2(
  RVG5_venom_circos_df,
  max_value = max_miRNA_count,
  track_height = trackheight,
  chromosome_col = 'miRNA.Sequence.Chrom',
  loci_location_col = 'miRNA.Start',
  value_col = 'miRNA.RPM.RVG.5S.CV1087.viridis.North.F',
  color_ramp = miRNA_color_ramp,
  chrom_labels = F
)

# Call the function to draw lines
draw_target_arrows(
  RVG5_venom_circos_df,
  sector1_col = 'miRNA.Sequence.Chrom',
  position1_col = 'miRNA.Start',
  sector2_col = 'miRNA.Target.Chrom',
  position2_col = 'miRNA.Target.Position',
  line_color_col = 'Color',
  line_type_col = 'Line.Type'
)

# Finish pdf
dev.off()


#### CIRCOS Plot for RVG_6S colored by venom gene family ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/CIRCOS_Functions.R')


# Create smaller data frame for circos plot circos.link function to draw lines
# This data frame also only has RVG_6S information in it and is color coded for venom gene families
RVG6_venom_circos_df <- venom_circos_df %>% 
  dplyr::select(
    -contains('.RVG.5S.'),
    -contains('.RVG.7S.'),
    -contains('.LVG.2.'),
    -contains('.LVG.4.'),
    -contains('.LVG.9.')
  ) %>% 
  filter(!miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M == 0) %>% # filter out any miRNA with zero expression, I don't think there are any, but just to be safe
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  distinct()

# Initialize max count values
max_protein_value <- range(log(RVG6_venom_circos_df$Intensity.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
max_miRNA_count <- range(log(RVG6_venom_circos_df$miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG6_venom_circos_df$RNA.VST.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)

# Define breaks for color scales
protein_breaks <- quantile(log(RVG6_venom_circos_df$Intensity.RVG.6S.CV0987.lutosus.Other.M + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
miRNA_breaks <- quantile(log(RVG6_venom_circos_df$miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
mRNA_breaks <- quantile(log(RVG6_venom_circos_df$RNA.VST.RVG.6S.CV0987.lutosus.Other.M + 1), probs = c(0, 0.5, 1), na.rm = TRUE)

# Create color ramps using magma color map
protein_color_ramp <- colorRamp2(protein_breaks, magma(3))
miRNA_color_ramp <- colorRamp2(miRNA_breaks, magma(3))
mRNA_color_ramp <- colorRamp2(mRNA_breaks, magma(3))


# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/Filtered/Venom_Genes_Colored_CIRCOS_CV0987_lutosus_2024.09.02.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par('track.height' = trackheight,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = c(rep(1, 17), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# # Draw circos tracks for protein
# draw_circos_tracks_points2(
#   RVG6_venom_circos_df,
#   max_value = max_protein_value,
#   track_height = trackheight,
#   chromosome_col = 'miRNA.Target.Chrom',
#   loci_location_col =  'miRNA.Target.Start',
#   value_col = 'Intensity.RVG.6S.CV0987.lutosus.Other.M',
#   color_ramp = protein_color_ramp,
#   chrom_labels = T
#   )

# Draw circos tracks for mRNA
draw_circos_tracks_points2(
  RVG6_venom_circos_df,
  max_value = max_mRNA_expression,
  track_height = trackheight,
  chromosome_col = 'miRNA.Target.Chrom',
  loci_location_col = 'miRNA.Target.Start',
  value_col = 'RNA.VST.RVG.6S.CV0987.lutosus.Other.M',
  color_ramp = mRNA_color_ramp,
  chrom_labels = T
)

# Draw circos tracks for miRNA
draw_circos_tracks_points2(
  RVG6_venom_circos_df,
  max_value = max_miRNA_count,
  track_height = trackheight,
  chromosome_col = 'miRNA.Sequence.Chrom',
  loci_location_col = 'miRNA.Start',
  value_col = 'miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M',
  color_ramp = miRNA_color_ramp,
  chrom_labels = F
)


# Call the function to draw lines
draw_target_arrows(RVG6_venom_circos_df,
                   sector1_col = 'miRNA.Sequence.Chrom',
                   position1_col = 'miRNA.Start',
                   sector2_col = 'miRNA.Target.Chrom',
                   position2_col = 'miRNA.Target.Position',
                   line_color_col = 'Color',
                   line_type_col = 'Line.Type'
                   )

# Finish pdf
dev.off()




#### Mini CICROS plot for PLA2 genes only ####
# 
# 
# # Load CIRCOS plot packages
# library(circlize)
# library(rlist)
# library(officer)
# library(rvg)
# source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/CIRCOS_Functions.R')
# 
# # Create smaller data frame for circos plot circos.link function to draw lines
# PLA2_circos_df <- mi_df %>% 
#   dplyr::select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, Positions,
#          miRNA.Counts.RVG.5S.CV1087.viridis.North.F, miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M, miRNA.Counts.RVG.7S.CV0985.concolor.Other.F,
#          RNA.VST.RVG.5S.CV1087.viridis.North.F, RNA.VST.RVG.6S.CV0987.lutosus.Other.M, RNA.VST.RVG.7S.CV0985.concolor.Other.F,
#          Intensity.RVG.5S.CV1087.viridis.North.F, Intensity.RVG.6S.CV0987.lutosus.Other.M, Intensity.RVG.7S.CV0985.concolor.Other.F) %>%
#   filter(str_detect(Genes, 'PLA2')) %>% # Filter out everything but PLA2
#   #arrange(desc(Intensity.RVG.5S.CV1087.viridis.North.F)) %>%
#   filter(!miRNA.Target.Chrom %in% c('PE-reconstructed-10x-myo', 'scaffold-un187')) %>%
#   filter(!miRNA.Sequence.Chrom %in% c('scaffold-un11', 'scaffold-un619', 'scaffold-un31', 'scaffold-un147')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
#   mutate(Color = case_when(
#     str_starts(Genes, 'Venom_PLA2') ~ 'purple',
#     str_starts(Genes, 'Venom_SVSP') ~ 'blue',
#     str_starts(Genes, 'Venom_SVMP') ~ 'red',
#     str_starts(Genes, 'Venom_vQC') ~ 'violet',
#     str_starts(Genes, 'Venom_ohanin') ~ 'green',
#     str_starts(Genes, 'Venom_LAAO') ~ 'orange',
#     str_starts(Genes, 'Venom_CTL') ~ 'magenta',
#     str_starts(Genes, 'Venom_CRISP') ~ 'coral',
#     str_starts(Genes, 'Venom_BPP') ~ 'cyan',
#     str_starts(Genes, 'Venom_VEGF') ~ 'navy',
#     str_starts(Genes, 'Venom_myotoxin') ~ 'maroon',
#     str_starts(Genes, 'Venom_BPP') ~ 'gray',
#     str_starts(Genes, 'Venom_EXO') ~ 'pink',
#     str_starts(Genes, 'Venom_ADAM') ~ 'aquamarine',
#     T ~ 'black'
#   ))
# 
# 
# # Create individualized circos data frames for each sample
# RVG5_PLA2_circos_df <- PLA2_circos_df %>% 
#   dplyr::select(-miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M, -miRNA.Counts.RVG.7S.CV0985.concolor.Other.F,
#          -RNA.VST.RVG.6S.CV0987.lutosus.Other.M, -RNA.VST.RVG.7S.CV0985.concolor.Other.F,
#          -Intensity.RVG.6S.CV0987.lutosus.Other.M, -Intensity.RVG.7S.CV0985.concolor.Other.F
#   ) %>% 
#   dplyr::rename(
#     miRNA = miRNA.Counts.RVG.5S.CV1087.viridis.North.F,
#     Protein = Intensity.RVG.5S.CV1087.viridis.North.F,
#     RNA = RNA.VST.RVG.5S.CV1087.viridis.North.F
#   ) %>%
#   mutate(miRNA.Target.End = miRNA.Target.End - miRNA.Target.Start) %>%
#   mutate(miRNA.Target.Start = 0) %>% 
#   mutate(miRNA.End = miRNA.End - miRNA.Start) %>% 
#   mutate(miRNA.Start = 0) %>% 
#   #dplyr::select(-miRNA.Cluster, -miRNA.Sequence.Chrom, -miRNA) %>% 
#   distinct()
# 
# # Reset the row names of the above dataframe
# rownames(RVG5_PLA2_circos_df) <- NULL
# 
# # Remove and alter only the miRNA loci from the dataframe
# miRNA_scaffold_df <- RVG5_PLA2_circos_df %>% 
#   dplyr::select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End) %>% 
#   dplyr::rename(
#     Loci = miRNA.Cluster,
#     Chrom = miRNA.Sequence.Chrom,
#     Start = miRNA.Start,
#     End = miRNA.End
#   )
# 
# # Now take only the target loci information
# target_scaffold_df <- RVG5_PLA2_circos_df %>% 
#   dplyr::select(Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End) %>% 
#   dplyr::rename(
#     Loci = Genes,
#     Chrom = miRNA.Target.Chrom,
#     Start = miRNA.Target.Start,
#     End = miRNA.Target.End
#   )
# 
# # Fuse the bottom of the two above data frames so that I now have a scaffold file I can use to draw sectors
# scaffold_pla2_df <- bind_rows(miRNA_scaffold_df, target_scaffold_df) %>% 
#   dplyr::rename(Size = End) %>% 
#   distinct()
# 
# # Reset indexing because it all weird and now:
# rownames(scaffold_pla2_df) <- NULL
# 
# # Initialize max count values
# max_protein_value <- range(log(RVG5_PLA2_circos_df$Protein + 1), na.rm = T)
# max_miRNA_count <- range(log(RVG5_PLA2_circos_df$miRNA + 1), na.rm = T)
# max_mRNA_expression <- range(log(RVG5_PLA2_circos_df$RNA + 1), na.rm = T)
# 
# # Create filtere
# 
# # Set sector height:
# sector_height <- 0.2
# 
# # Begin saving the plot as a pdf
# pdf(file = 'Figures/CIRCOS/3UTR/PLA_mini_CIRCOS.2024.09.02.pdf', width = 10, height = 10)
# 
# circos.clear()
# 
# # Set basic parameters for the circos plot
# circos.par(track.height = sector_height,
#            cell.padding = c(0.01, 0, 0.01, 0),
#            start.degree = 90,
#            gap.degree = 15)
# #gap.degree = c(rep(1, 15), 20))
# 
# # Initiate the circos plot
# # xlim is used to specify the start and the end value of each chromosome (sectors)
# circos.initialize(scaffold_pla2_df$Loci, xlim = cbind(scaffold_pla2_df$Start, scaffold_pla2_df$Size))
# 
# # Create data frame with only miRNA information
# mi_expression_df <- RVG5_PLA2_circos_df %>% 
#   dplyr::select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, miRNA) %>% 
#   distinct() %>% 
#   dplyr::rename(
#     Loci = miRNA.Cluster,
#     Chrom = miRNA.Sequence.Chrom,
#     Start = miRNA.Start,
#     End = miRNA.End,
#     Expression = miRNA
#   ) %>% 
#   mutate(
#     Normalized.Expression = scale(Expression)
#   )
# 
# # Now for a data frame for only geneic information
# gene_expression_df <- RVG5_PLA2_circos_df %>% 
#   dplyr::select(Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, RNA) %>% 
#   distinct() %>% 
#   dplyr::rename(
#     Loci = Genes,
#     Chrom = miRNA.Target.Chrom,
#     Start = miRNA.Target.Start,
#     End = miRNA.Target.End,
#     Expression = RNA
#   ) %>% 
#   mutate(
#     Normalized.Expression = scale(Expression)
#   )
# 
# 
# # Define new circos data frame that combines all of the expression into a single column
# combined_expression_df <- bind_rows(mi_expression_df, gene_expression_df) %>% 
#   distinct()
# 
# circos.track(
#   ylim = c(0, sector_height),
#   track.height = sector_height,
#   panel.fun = function(x, y) {
#     
#     # Filter data for current chromosome
#     # Create variable that holds the sector data using the genes as sectors
#     sector <- combined_expression_df$Loci
#     
#     # Create variable to hold the information for chromomes
#     chrom <- combined_expression_df$Chrom
#     
#     # Variable for the current sector
#     current_sector <- combined_expression_df %>% 
#       filter(sector == CELL_META$sector.index)
#     
#     # Iterate over the rows of current_sector
#     for (loci in seq_len(nrow(current_sector))) {
#       
#       # Variables for x and y
#       start <- current_sector$Start[loci]
#       end <- current_sector$End[loci]
#       value <- log(current_sector$Expression[loci] + 1)
#       norm_value <- current_sector$Normalized.Expression[loci]
#       
#       # # Assign color based on value
#       # color <- heat.colors(10)[findInterval(value, seq(min(max_mRNA_expression), max(max_mRNA_expression), length.out = 10))]
#       # Assign color based on normalized value
#       color <- heat.colors(10)[findInterval(norm_value, seq(-3, 3, length.out = 10))]
#       
#       # Debugging
#       print(start)
#       print(end)
#       print(value)
#       print(norm_value)
#       
#       # Draw rectangles
#       circos.rect(xleft = start,
#                   xright = end,
#                   ybottom = 0,
#                   ytop = sector_height,
#                   col = color)
#       
#       # Add text to the CIRCOS that represents the chromosome and the gene
#       circos.text(
#         CELL_META$xcenter,
#         CELL_META$cell.ylim[2] + mm_y(7),
#         gsub('Venom_', '', CELL_META$sector.index),
#         facing = 'downward',
#         cex = 0.6
#       )
#     }
#   }
# )
# 
# # Again create a new data frame for the data so that I can draw arrows
# arrows_PLA2_df <- RVG5_PLA2_circos_df %>% 
#   mutate(miRNA.Midpoint = (miRNA.Start + miRNA.End) / 2) %>% 
#   mutate(Positions = as.numeric(Positions))
# 
# # Plot each miRNA target line one by one with a for loop
# for(cluster in 1:nrow(arrows_PLA2_df)) {
#   
#   # Need to specify sector1, position1 and sector2, position2 to plot the miRNA targeting line
#   sector1 <- arrows_PLA2_df$miRNA.Cluster[cluster]
#   position1 <- arrows_PLA2_df$miRNA.Midpoint[cluster]
#   sector2 <- arrows_PLA2_df$Genes[cluster]
#   position2 <- arrows_PLA2_df$Positions[cluster]
#   line_and_arrow_color <- arrows_PLA2_df$Color[cluster]
#   
#   # Print for debugging
#   print(arrows_PLA2_df[cluster,])
#   
#   # Draw links between sectors based on sector and position.
#   circos.link(
#     sector1,
#     position1,
#     sector2,
#     position2,
#     col = line_and_arrow_color,
#     #col = add_transparency(line_and_arrow_color, 0.95), # Adjust the line color as needed
#     directional = 1,
#     lwd = 1,
#     arr.width = 0.1,
#     arr.length = 0.1
#   )
# }
# 
# # Save
# dev.off()



#### CIRCOS Plot for RVG_7S colored by venom gene family ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/CIRCOS_Functions.R')


# Create smaller data frame for circos plot circos.link function to draw lines
# This data frame also only has RVG_7S information in it and is color coded for venom gene families
RVG7_venom_circos_df <- venom_circos_df %>% 
  dplyr::select(
    -contains('.RVG.5S.'),
    -contains('.RVG.6S.'),
    -contains('.LVG.2.'),
    -contains('.LVG.4.'),
    -contains('.LVG.9.')
  ) %>% 
  filter(!miRNA.RPM.RVG.7S.CV0985.concolor.Other.F == 0) %>% # filter out any miRNA with zero expression, I don't think there are any, but just to be safe
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  distinct()

# Initialize max count values
max_protein_value <- range(log(RVG7_venom_circos_df$Intensity.RVG.7S.CV0985.concolor.Other.F + 1), na.rm = T)
max_miRNA_count <- range(log(RVG7_venom_circos_df$miRNA.RPM.RVG.7S.CV0985.concolor.Other.F + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG7_venom_circos_df$RNA.VST.RVG.7S.CV0985.concolor.Other.F + 1), na.rm = T)

# Define breaks for color scales
protein_breaks <- quantile(log(RVG7_venom_circos_df$Intensity.RVG.7S.CV0985.concolor.Other.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
miRNA_breaks <- quantile(log(RVG7_venom_circos_df$miRNA.RPM.RVG.7S.CV0985.concolor.Other.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
mRNA_breaks <- quantile(log(RVG7_venom_circos_df$RNA.VST.RVG.7S.CV0985.concolor.Other.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)

# Create color ramps using magma color map
protein_color_ramp <- colorRamp2(protein_breaks, magma(3))
miRNA_color_ramp <- colorRamp2(miRNA_breaks, magma(3))
mRNA_color_ramp <- colorRamp2(mRNA_breaks, magma(3))

# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/Filtered/Venom_Genes_Colored_CIRCOS_CV0985_concolor_2024.09.02.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par('track.height' = trackheight,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = c(rep(1, 17), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# # Draw circos tracks for protein
# draw_circos_tracks_points2(
#   RVG7_venom_circos_df,
#   max_value = max_protein_value,
#   track_height = trackheight,
#   chromosome_col = 'miRNA.Target.Chrom',
#   loci_location_col =  'miRNA.Target.Start',
#   value_col = 'Intensity.RVG.7S.CV0985.concolor.Other.F',
#   color_ramp = protein_color_ramp,
#   chrom_labels = T
#   )

# Draw circos tracks for mRNA
draw_circos_tracks_points2(
  RVG7_venom_circos_df,
  max_value = max_mRNA_expression,
  track_height = trackheight,
  chromosome_col = 'miRNA.Target.Chrom',
  loci_location_col = 'miRNA.Target.Start',
  value_col = 'RNA.VST.RVG.7S.CV0985.concolor.Other.F',
  color_ramp = mRNA_color_ramp,
  chrom_labels = T
)

# Draw circos tracks for miRNA
draw_circos_tracks_points2(
  RVG7_venom_circos_df,
  max_value = max_miRNA_count,
  track_height = trackheight,
  chromosome_col = 'miRNA.Sequence.Chrom',
  loci_location_col = 'miRNA.Start',
  value_col = 'miRNA.RPM.RVG.7S.CV0985.concolor.Other.F',
  color_ramp = miRNA_color_ramp,
  chrom_labels = F
)


# Call the function to draw lines
draw_target_arrows(RVG7_venom_circos_df,
                   sector1_col = 'miRNA.Sequence.Chrom',
                   position1_col = 'miRNA.Start',
                   sector2_col = 'miRNA.Target.Chrom',
                   position2_col = 'miRNA.Target.Position',
                   line_color_col = 'Color',
                   line_type_col = 'Line.Type'
                   )

# Finish pdf
dev.off()





#### CIRCOS All Samples ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)


# Rename RNA.VST to mRNA.VST
mi_df <- mi_df %>% 
  rename_with(
    ~ sub('^RNA.VST', 'mRNA.VST', .), starts_with('RNA.VST')
  )


# Create smaller data frame for circos plot circos.link function to draw lines
venom_circos_df2 <- mi_df %>% 
  dplyr::select(
    miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, Positions, Origin,
    -contains('VST'),
    -contains('Intensity.'),
    -contains('RPM'),
    -contains('Counts')
  ) %>% 
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  filter(!miRNA.Target.Chrom %in% c('PE-reconstructed-10x-myo')) %>% 
  filter(!str_detect(miRNA.Sequence.Chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  filter(!str_detect(miRNA.Target.Chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  mutate(Positions = strsplit(as.character(Positions), ' ')) %>% # Separate the Positions column into multiple rows
  unnest(Positions) %>% 
  mutate(Positions = as.numeric(Positions)) %>%  # Convert Positions to numeric
  filter(!is.na(Positions)) %>% 
  # Create the new column miRNA.Target.Position
  mutate(miRNA.Target.Position = miRNA.Target.Start + Positions) %>% 
  mutate(Venom.Color = case_when(
    str_starts(Genes, 'Venom_PLA2') ~ PLA2_color,
    str_starts(Genes, 'Venom_SVSP') ~ SVSP_color,
    str_starts(Genes, 'Venom_SVMP') ~ SVMP_color,
    str_starts(Genes, 'Venom_ADAM') ~ ADAM_color,
    str_starts(Genes, 'Venom_vQC') ~ vQC_color,
    str_starts(Genes, 'Venom_ohanin') ~ ohanin_color,
    str_starts(Genes, 'Venom_LAAO') ~ LAAO_color,
    str_starts(Genes, 'Venom_CTL') ~ CTL_color,
    str_starts(Genes, 'Venom_CRISP') ~ CRISP_color,
    str_starts(Genes, 'Venom_BPP') ~ BPP_color,
    str_starts(Genes, 'Venom_VEGF') ~ VEGF_color,
    str_starts(Genes, 'Venom_myotoxin') ~ myotoxin_color,
    str_starts(Genes, 'Venom_EXO') ~ EXO_color,
    T ~ other_color
  )) %>%
  mutate(
    miRNA.Color = miRNA_color
  ) %>% 
  mutate(
    Target.Color = case_when(
      str_detect(Origin, 'three_prime_utr') ~ three_prime_color,
      str_detect(Origin, 'five_prime_utr') ~ five_prime_color,
      str_detect(Origin, 'CDS') ~ cds_color
    )
  ) %>% 
  distinct()



# Create another data frame that only contains miRNA data for later fusion
mirna_df <- mi_df %>% 
  dplyr::select(
    miRNA.Cluster, contains('miRNA.VST')
  ) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('miRNA.VST'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.VST'
  ) %>% 
  group_by(miRNA.Cluster) %>% 
  summarize(
    miRNA.VST.Mean = mean(miRNA.VST),
    miRNA.VST.Variance = var(miRNA.VST)
  )

# Create another data frame that only contains mRNA data for later fusion
mrna_df <- mi_df %>% 
  dplyr::select(
    Genes, contains('mRNA.VST')
  ) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('mRNA.VST'),
    names_to = 'Sample.ID',
    values_to = 'mRNA.VST'
  ) %>%
  group_by(Genes) %>% 
  summarize(
    mRNA.VST.Mean = mean(mRNA.VST),
    mRNA.VST.Variance = var(mRNA.VST)
  )

# Create another data frame that only contains protein data for later fusion
protein_df <- mi_df %>% 
  dplyr::select(
    Genes, contains('Intensity')
  ) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('Intensity'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>% 
  group_by(Genes) %>% 
  summarize(
    Intensity.Mean = mean(Intensity),
    Intensity.Var = var(Intensity)
  )


# Fuse the data frames back together
venom_circos_df3 <- left_join(
  venom_circos_df2,
  mirna_df,
  by = c('miRNA.Cluster')
) %>% 
  left_join(
    mrna_df,
    by = c('Genes')
  ) %>% 
  left_join(
    protein_df,
    by = c('Genes')
  )




# START THE DRAWING Proccess

# Initialize max count values
max_protein_value <- range(log(venom_circos_df3$Intensity.Mean + 1), na.rm = T)
max_protein_variance <- range(log(venom_circos_df3$Intensity.Mean + 1) + sqrt(log(venom_circos_df3$Intensity.Var + 1)), na.rm = T)
ylim_protein <- c(min(max_protein_value), max(max_protein_variance))

max_mRNA_expression <- range(venom_circos_df3$mRNA.VST.Mean, na.rm = T)
max_mRNA_variance <- range(venom_circos_df3$mRNA.VST.Mean + sqrt(venom_circos_df3$mRNA.VST.Variance), na.rm = T)
ylim_mRNA <- c(min(max_mRNA_expression), max(max_mRNA_variance))

max_miRNA_expression <- range(venom_circos_df3$miRNA.VST.Mean, na.rm = T)
max_miRNA_variance <- range(venom_circos_df3$miRNA.VST.Mean + sqrt(venom_circos_df3$miRNA.VST.Variance), na.rm = T)
ylim_miRNA <- c(min(max_miRNA_expression), max(max_miRNA_variance))

# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)


# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/Filtered/All_Samples_CIRCOS_2024.09.11.pdf', width = 10, height = 10)

# Clear plot of any current CIRCOS plot
circos.clear()

# Set basic parameters for the circos plot
circos.par(
  track.height = trackheight,
  cell.padding = c(0.01, 0, 0.01, 0),
  start.degree = 90,
  gap.degree = c(rep(1, 17), 20)
)

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# Set track height
trackheight <- 0.2

# Set point size
cex <- 0.8

# Draw the protein expression points with variance bars
circos.track(
  ylim = ylim_protein,
  track.height = trackheight,
  panel.fun = function(x, y) {
    
    # Draw the protein points
    # Filter the data for the current chromosome
    current_chrom <- venom_circos_df3[venom_circos_df3['miRNA.Target.Chrom'] == CELL_META$sector.index, ]
    
    # Use foreach for iteration over rows
    foreach(loci = 1:nrow(current_chrom)) %do% {
      
      # Variables for x and y in the graph
      horizontal <- current_chrom[['miRNA.Target.Start']][loci]
      vertical <- log(current_chrom[['Intensity.Mean']][loci] + 1)
      variance <- log(current_chrom[['Intensity.Var']][loci] + 1)
      point_color <- current_chrom[['Venom.Color']][loci]
      
      # Calculate the upper and lower bounds for the confidence intervals
      lower_bound <- vertical - sqrt(variance)
      upper_bound <- vertical + sqrt(variance)
      
      # Print for debugging
      print(horizontal)
      print(vertical)
      
      # Draw points
      circos.points(
        horizontal,
        vertical,
        col = point_color,
        pch = 17,
        cex = cex
      )
      
      # Draw confidence intervals as error bars
      circos.segments(
        x0 = horizontal,
        x1 = horizontal,
        y0 = lower_bound,
        y1 = upper_bound,
        col = point_color
      )
    }
    
    # Add text to CIRCOS about chromosomes
    circos.text(
      CELL_META$xcenter,
      CELL_META$cell.ylim[2] + mm_y(7),
      gsub('scaffold-', '', CELL_META$sector.index),
      facing = 'downward',
      cex = cex
    )
  }
)


# Draw the mRNA and miRNA expression points with variance bars
circos.track(
  ylim = ylim_miRNA,
  track.height = trackheight,
  panel.fun = function(x, y) {
    
    # Draw the mRNA points
    # Filter the data for the current chromosome
    current_chrom1 <- venom_circos_df3[venom_circos_df3['miRNA.Target.Chrom'] == CELL_META$sector.index, ]
    
    # Use foreach for iteration over rows
    foreach(loci = 1:nrow(current_chrom1)) %do% {
      
      # Variables for x and y in the graph
      horizontal1 <- current_chrom1[['miRNA.Target.Start']][loci]
      vertical1 <- current_chrom1[['mRNA.VST.Mean']][loci]
      variance1 <- current_chrom1[['mRNA.VST.Variance']][loci]
      point_color1 <- current_chrom1[['Venom.Color']][loci]
      
      # Calculate the upper and lower bounds for the confidence intervals
      lower_bound1 <- vertical1 - sqrt(variance1)
      upper_bound1 <- vertical1 + sqrt(variance1)
      
      # Print for debugging
      print(horizontal1)
      print(vertical1)
      
      # Draw points
      circos.points(
        horizontal1,
        vertical1,
        col = point_color1,
        pch = 16,
        cex = cex
      )
      
      # Draw confidence intervals as error bars
      circos.segments(
        x0 = horizontal1,
        x1 = horizontal1,
        y0 = lower_bound1,
        y1 = upper_bound1,
        col = point_color1
      )
    }
    
    # Draw the miRNA points
    # Filter the data for the current chromosome
    current_chrom2 <- venom_circos_df3[venom_circos_df3['miRNA.Sequence.Chrom'] == CELL_META$sector.index, ]
    
    # Use foreach for iteration over rows
    foreach(loci = 1:nrow(current_chrom2)) %do% {
      
      # Variables for x and y in the graph
      horizontal2 <- current_chrom2[['miRNA.Start']][loci]
      vertical2 <- current_chrom2[['miRNA.VST.Mean']][loci]
      variance2 <- current_chrom2[['miRNA.VST.Variance']][loci]
      point_color2 <- current_chrom2[['miRNA.Color']][loci]
      
      # Calculate the upper and lower bounds for the confidence intervals
      lower_bound2 <- vertical2 - sqrt(variance2)
      upper_bound2 <- vertical2 + sqrt(variance2)
      
      # Print for debugging
      print(horizontal2)
      print(vertical2)
      
      # Draw points
      circos.points(
        horizontal2,
        vertical2,
        col = point_color2,
        pch = 18,
        cex = cex
      )
      
      # Draw confidence intervals as error bars
      circos.segments(
        x0 = horizontal2,
        x1 = horizontal2,
        y0 = lower_bound2,
        y1 = upper_bound2,
        col = point_color2
      )
    }
    
    # # Add text to CIRCOS about chromosomes
    # circos.text(
    #   CELL_META$xcenter,
    #   CELL_META$cell.ylim[2] + mm_y(7),
    #   gsub('scaffold-', '', CELL_META$sector.index),
    #   facing = 'downward',
    #   cex = 0.8
    # )
  }
)


# Create links between miRNA loci and mRNA Venom Genes
for (cluster in 1:nrow(venom_circos_df3)) {
  
  # Specify sectors, positions, and so on for miRNA targeting lines
  sector1 <- venom_circos_df3[['miRNA.Sequence.Chrom']][cluster] # Sector for the miRNAs
  position1 <- venom_circos_df3[['miRNA.Start']][cluster]
  sector2 <- venom_circos_df3[['miRNA.Target.Chrom']][cluster]
  position2 <- venom_circos_df3[['miRNA.Target.Start']][cluster]
  line_and_arrow_color <- venom_circos_df3[['Target.Color']][cluster]
  
  # Print for debugging
  print(venom_circos_df3[cluster, ])
  
  # Draw links between sectors based on sector and position
  circos.link(
    sector1,
    position1,
    sector2,
    position2,
    # col = line_and_arrow_color,
    col = add_transparency(line_and_arrow_color, 0.7), # Adjust the line color as needed
    directional = 1,
    lwd = 0.7,
    arr.width = 0.1,
    arr.length = 0.1
  )
}


# Assuming you have color schemes for miRNA clusters and Venom Gene families
miRNA_colors <- c("microRNA" = miRNA_color)
venom_gene_colors <- c(
  'SVMP' = SVMP_color,
  'ADAM' = ADAM_color,
  'SVSP' = SVSP_color, 
  'PLA2' = PLA2_color,
  'VEGF' = VEGF_color,
  'Ohanin' = ohanin_color,
  'Myotoxin' = myotoxin_color,
  'vQC' = vQC_color,
  'CRISP' = CRISP_color,
  'CTL' = CTL_color,
  'EXO' = EXO_color,
  'LAAO' = LAAO_color,
  'BPP' = BPP_color
)

# Combine both for the legend
all_colors <- c(miRNA_colors, venom_gene_colors)
labels <- names(all_colors)

# Add the legend
legend("topright", legend = labels, fill = all_colors, title = "Gene Family", cex = 0.8, bty = "n")

# Finish saving pdf
dev.off()
