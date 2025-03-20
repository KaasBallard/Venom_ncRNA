# Last Edited: 2025/02/05

# Set up and Read Data ----

## Load in packages ----
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
library(arrow)

## Load CIRCOS plot packages ----
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Functions/CIRCOS_Functions.R')


## Read in data ----
# Create variable for the fused data set.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read in the data
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)


## Math to calculate rough location of BPP in viridis based on adamanteus ----
# Length ma5 in both species and the location of bbp in C. adamanteus
bpp_adamanteus_start <- 101936827
bpp_adamanteus_end <- 101950398
adamanteus_ma5_length <- 106188508 # I got this from the .fai that samtools faidx GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_genomic.fna created. Note that the chromosome is CM077922.1
viridis_ma5_length <- 87326383

# Calculate the start position of bpp in viridis
bpp_viridis_start <- (viridis_ma5_length / adamanteus_ma5_length) * bpp_adamanteus_start


## Math to caculate the start site for myotoxin
nbisL4_mrna_4_start <- 289306016 # Bdef3.1
nbisL5_mrna_5_end <- 289328253 # Bdef4

# Find the midpoint of the two above bdef genes
myotoxin_start <- (nbisL4_mrna_4_start + nbisL5_mrna_5_end) / 2

## Format the data ----
# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    !is.na(feature.type)
  ) %>% # This should get rid of all of the weirdly annotated genes
  mutate(
    miRNA.target.length = miRNA.target.end - miRNA.target.start + 1
  ) %>% 
  # Add BPP start based on it's location in C. adamanteus
  mutate(
    miRNA.target.start = case_when(
      str_detect(genes, 'Venom_BPP') ~ bpp_viridis_start,
      str_detect(genes, 'Venom_myotoxin') ~ myotoxin_start,
      TRUE ~ miRNA.target.start
    )
  ) %>% 
  # Get the BPP end based off of length and start
  mutate(
    miRNA.target.end = case_when(
      str_detect(genes, 'Venom_BPP') ~ miRNA.target.start + miRNA.target.length,
      str_detect(genes, 'Venom_myotoxin') ~ miRNA.target.start + miRNA.target.length,
      TRUE ~ miRNA.target.end
    )
  ) %>% 
  # Correct the chrom for myotoxin and BPP
  mutate(
    miRNA.target.chrom = case_when(
      str_detect(genes, 'Venom_BPP') ~ 'scaffold-ma5', # Add correct BPP chrom
      str_detect(genes, 'Venom_myotoxin') ~ 'scaffold-ma1', # Add correct myotoxin chrom
      TRUE ~ miRNA.target.chrom
    )
  ) 
rm(miRNA_mRNA_protein_df)


# Add more filters the main data frame separately so things don't get to complicated
mi_df <- mi_df %>% 
  filter(total.energy <= -7) %>% # I am including this filter to remove any low strength binding examples
  filter(total.score >= 155) %>% # I am including this filter to remove any low confidence binding examples
  mutate(
    line.type = case_when(
      feature.type == 'three_prime_utr' ~ 'solid',
      feature.type == 'five_prime_utr' ~ 'longdash',
      feature.type == 'CDS' ~ 'twodash',
    )
  )
  

# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
SVSP_color <- '#F0B830' 
PLA2_color <- '#7570B3'
VEGF_color <- '#74ADD1'
ohanin_color <- '#3A489C'
myotoxin_color <- '#B2182B'
vQC_color <- '#80BC50'
CRISP_color <- '#E7298A'
CTL_color <- '#F67E17'
EXO_color <- '#005824'
LAAO_color <- '#B35806'
BPP_color <- '#1B9E77'
other_color <- '#666666'
miRNA_color <- 'grey'
three_prime_color <- 'black'
# five_prime_color <- '#0072b2'
five_prime_color <- 'red'
# cds_color <- '#d55e00'
# cds_color <- '#4A70B5'
cds_color <- 'grey'


# Create General Dataframe for CIRCOS ----

# Create smaller data frame for circos plot circos.link function to draw lines
venom_circos_df <- mi_df %>% 
  dplyr::select(
    miRNA.cluster, miRNA.sequence.chrom, miRNA.start, miRNA.end, genes, miRNA.target.chrom, miRNA.target.start, miRNA.target.end, positions, feature.type, line.type
  ) %>% 
  filter(str_starts(genes, 'Venom'), !str_detect(genes, 'ADAM')) %>% # Filter out everything but venom
  filter(!miRNA.target.chrom %in% c('PE-reconstructed-10x-myo')) %>% 
  filter(!str_detect(miRNA.sequence.chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  filter(!str_detect(miRNA.target.chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  mutate(positions = strsplit(as.character(positions), ' ')) %>% # Separate the positions column into multiple rows
  unnest(positions) %>% 
  mutate(positions = as.numeric(positions)) %>% # Convert positions to numeric
  filter(!is.na(positions)) %>% 
  mutate(
    # Create the new column miRNA.target.position
    miRNA.target.position = miRNA.target.start + positions,
    # Create colors for the venoms
    venom.color = case_when(
      str_starts(genes, 'Venom_PLA2') ~ PLA2_color,
      str_starts(genes, 'Venom_SVSP') ~ SVSP_color,
      str_starts(genes, 'Venom_SVMP') ~ SVMP_color,
      str_starts(genes, 'Venom_vQC') ~ vQC_color,
      str_starts(genes, 'Venom_ohanin') ~ ohanin_color,
      str_starts(genes, 'Venom_LAAO') ~ LAAO_color,
      str_starts(genes, 'Venom_CTL') ~ CTL_color,
      str_starts(genes, 'Venom_CRISP') ~ CRISP_color,
      str_starts(genes, 'Venom_BPP') ~ BPP_color,
      str_starts(genes, 'Venom_VEGF') ~ VEGF_color,
      str_starts(genes, 'Venom_myotoxin') ~ myotoxin_color,
      str_starts(genes, 'Venom_EXO') ~ EXO_color,
      T ~ other_color
    ),
    # miRNA color
    miRNA.color = miRNA_color,
    # Create colors for the targets types
    target.color = case_when(
      str_detect(feature.type, 'three_prime_utr') ~ three_prime_color,
      str_detect(feature.type, 'five_prime_utr') ~ five_prime_color,
      str_detect(feature.type, 'CDS') ~ cds_color
    )
  ) %>%
  distinct()

# Create another data frame that only contains miRNA data for later fusion
mirna_df <- mi_df %>% 
  dplyr::select(
    sample.id, miRNA.cluster, miRNA.vst
  ) %>% 
  distinct() %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    miRNA.vst.mean = mean(miRNA.vst),
    miRNA.vst.var = var(miRNA.vst)
  )

# Create another data frame that only contains mRNA data for later fusion
mrna_df <- mi_df %>% 
  dplyr::select(
    sample.id, genes, mRNA.vst
  ) %>% 
  distinct() %>% 
  group_by(genes) %>% 
  summarise(
    mRNA.vst.mean = mean(mRNA.vst),
    mRNA.vst.var = var(mRNA.vst)
  )

# Create another data frame that only contains protein data for later fusion
protein_df <- mi_df %>% 
  dplyr::select(
    sample.id, genes, intensity
  ) %>% 
  distinct() %>% 
  group_by(genes) %>% 
  summarise(
    intensityMean = mean(intensity),
    intensityVar = var(intensity)
  )


# Fuse the data frames back together
venom_circos_df2 <- left_join(
  venom_circos_df,
  mirna_df,
  by = c('miRNA.cluster')
) %>% 
  left_join(
    mrna_df,
    by = c('genes')
  ) %>% 
  left_join(
    protein_df,
    by = c('genes')
  )
glimpse(venom_circos_df2)


# CIRCOS All Samples ----

## Start drawing ----

# Initialize max count values
max_protein_value <- range(log(venom_circos_df2$intensityMean + 1), na.rm = T)
max_protein_variance <- range(log(venom_circos_df2$intensityMean + 1) + sqrt(log(venom_circos_df2$intensityVar + 1)), na.rm = T)
ylim_protein <- c(min(max_protein_value), max(max_protein_variance))

max_mRNA_expression <- range(venom_circos_df2$mRNA.vst.mean, na.rm = T)
max_mRNA_variance <- range(venom_circos_df2$mRNA.vst.mean + sqrt(venom_circos_df2$mRNA.vst.var), na.rm = T)
ylim_mRNA <- c(min(max_mRNA_expression), max(max_mRNA_variance))

max_miRNA_expression <- range(venom_circos_df2$miRNA.vst.mean, na.rm = T)
max_miRNA_variance <- range(venom_circos_df2$miRNA.vst.mean + sqrt(venom_circos_df2$miRNA.vst.var), na.rm = T)
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

# Set track height
trackheight <- 0.2

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/Filtered/All_Samples_CIRCOS_2025.02.05.pdf', width = 10, height = 10)

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

# Set point size
cex <- 0.8

# Draw the protein expression points with variance bars
circos.track(
  ylim = ylim_protein,
  track.height = trackheight,
  panel.fun = function(x, y) {
    
    # Draw the protein points
    # Filter the data for the current chromosome
    current_chrom <- venom_circos_df2[venom_circos_df2['miRNA.target.chrom'] == CELL_META$sector.index, ]
    
    # Use foreach for iteration over rows
    foreach(loci = 1:nrow(current_chrom)) %do% {
      
      # Variables for x and y in the graph
      horizontal <- current_chrom[['miRNA.target.start']][loci]
      vertical <- log(current_chrom[['intensityMean']][loci] + 1)
      variance <- log(current_chrom[['intensityVar']][loci] + 1)
      point_color <- current_chrom[['venom.color']][loci]
      
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
        pch = 16,
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
    current_chrom1 <- venom_circos_df2[venom_circos_df2['miRNA.target.chrom'] == CELL_META$sector.index, ]
    
    # Use foreach for iteration over rows
    foreach(loci = 1:nrow(current_chrom1)) %do% {
      
      # Variables for x and y in the graph
      horizontal1 <- current_chrom1[['miRNA.target.start']][loci]
      vertical1 <- current_chrom1[['mRNA.vst.mean']][loci]
      variance1 <- current_chrom1[['mRNA.vst.var']][loci]
      point_color1 <- current_chrom1[['venom.color']][loci]
      
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
    current_chrom2 <- venom_circos_df2[venom_circos_df2['miRNA.sequence.chrom'] == CELL_META$sector.index, ]
    
    # Use foreach for iteration over rows
    foreach(loci = 1:nrow(current_chrom2)) %do% {
      
      # Variables for x and y in the graph
      horizontal2 <- current_chrom2[['miRNA.start']][loci]
      vertical2 <- current_chrom2[['miRNA.vst.mean']][loci]
      variance2 <- current_chrom2[['miRNA.vst.var']][loci]
      point_color2 <- current_chrom2[['miRNA.color']][loci]
      
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
        pch = 16,
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


# Create links between miRNA loci and mRNA Venom genes
for (cluster in 1:nrow(venom_circos_df2)) {
  
  # Specify sectors, positions, and so on for miRNA targeting lines
  sector1 <- venom_circos_df2[['miRNA.sequence.chrom']][cluster] # Sector for the miRNAs
  position1 <- venom_circos_df2[['miRNA.start']][cluster]
  sector2 <- venom_circos_df2[['miRNA.target.chrom']][cluster]
  position2 <- venom_circos_df2[['miRNA.target.start']][cluster]
  line_and_arrow_color <- venom_circos_df2[['target.color']][cluster]
  
  # Print for debugging
  print(venom_circos_df2[cluster, ])
  
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
