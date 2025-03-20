# Last Edited: 2024/9/13

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
library(scales)

# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")

# Create variable for the fused dataset.
# If this is in Dropbox the file needs to be on the machine
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

# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
ADAM_color <- '#9A70B5'
SVSP_color <- '#F0B830' 
PLA2_color <- '#7570B3'
miRNA_color <- '#8B0AA5'
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
three_prime_color <- 'black'
# five_prime_color <- '#0072b2'
five_prime_color <- '#1B9E77'
# cds_color <- '#d55e00'
cds_color <- '#4A70B5'





#### Gene family plots dataframe for Venom Gene Dot Plots (VST) ####

# Create a data frame for venom genes only
venom_genes_vst_df <- mi_df %>%
  filter(
    str_starts(Genes, 'Venom_'),
  ) %>%
  filter(In.Library == 'Yes') %>% # Remove any proteins not in the library sent to Anthony
  dplyr::select(
    miRNA.Cluster,
    Genes,
    contains('miRNA.Counts.'),
    contains('RNA.VST.'),
    contains('Intensity.'),
    contains('miRNA.RPM.')
  ) %>%
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  grepl('VEGF', Genes) ~ 'VEGF',
                                  grepl('ohanin', Genes) ~ 'ohanin',
                                  grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('ADAM', Genes) ~ 'ADAM',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  grepl('CTL', Genes) ~ 'CTL',
                                  grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>%
  distinct()


#### mRNA (VST) vs Protein for all samples ####


# Pivot data frame so that all samples are in a single column
all_samples_plot_protein_df <- venom_genes_vst_df %>%
  dplyr::select(-contains('miRNA.'), -contains('Intensity')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('RNA.VST.'),
    names_to = 'Sample.ID',
    values_to = 'RNA.VST'
  ) %>%
  mutate(Sample.ID = gsub('RNA.VST.', '', Sample.ID))

# Pivot this one too
all_samples_plot_mRNA_df <- venom_genes_vst_df %>%
  dplyr::select(-contains('miRNA.'), -contains('RNA.')) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('Intensity.'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>% 
  mutate(Sample.ID = gsub('Intensity.', '', Sample.ID))

# Fuse back together
all_samples_plot_df <- left_join(
  all_samples_plot_protein_df,
  all_samples_plot_mRNA_df,
  by = c('Sample.ID', 'Genes', 'Venom.Family')
) 

# write.table(all_samples_plot_df, file = '/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Fxn_NSF/_Manuscripts/3_Venom_ncRNA/Data/Normalized_RNA_Protein.tsv', sep = '\t')


# Calculate CLR
all_samples_clr_df <- all_samples_plot_df %>%
  group_by(Sample.ID) %>%
  mutate(Scaled.Protein.Expression = Intensity / sum(Intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.mRNA.Expression = RNA.VST / sum(RNA.VST)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  mutate(CLR.Scaled.Protein.Expression = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(CLR.Scaled.mRNA.Expression = as.numeric(compositions::clr(Scaled.mRNA.Expression)))


# Set variables for protein and RNA expression levels
all_mRNA_expression <- all_samples_clr_df$CLR.Scaled.mRNA.Expression
all_protein_expression <- all_samples_clr_df$CLR.Scaled.Protein.Expression
venom_family <- all_samples_clr_df$Venom.Family
# gene_label <- all_samples_clr_df$Genes


# Create color scheme for the venom genes
venom_colors <- c(
  SVMP = '#4A70B5',
  ADAM = '#9A70B5',
  SVSP = '#F0B830', 
  PLA2 = '#7570B3',
  miRNA = '#8B0AA5',
  VEGF = '#74ADD1',
  ohanin = '#3A489C',
  myotoxin = '#B2182B',
  vQC = '#80BC50',
  CRISP = '#E7298A',
  CTL = '#F67E17',
  EXO = '#49FFFF',
  LAAO = '#B35806',
  BPP = '#1B9E77',
  others = '#666666'
)


# Find linear regression equation
venom_lm <- lm(all_protein_expression ~ all_mRNA_expression, data = all_samples_clr_df)
venom_slope <- coef(venom_lm)[2]
venom_intercept <- coef(venom_lm)[1]

# Create figure for mRNA and Protein expression
all_gene_expression_plot <- all_samples_clr_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
  geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +  # Apply color scheme  
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
all_gene_expression_plot
# It seems the data, at least if you plot it for all of the individuals is Heteroskedastic and violates the assumptions of the linear model

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = all_gene_expression_plot, create.dir = T)



#### Dataframe for residuals and miRNA (VST) ####

# Find linear regression equation
venom_lm <- lm(all_protein_expression ~ all_mRNA_expression, data = all_samples_clr_df)
venom_slope <- coef(venom_lm)[2]
venom_intercept <- coef(venom_lm)[1]

# Find the predicted values from the model and put them into the dataframe
all_samples_residuals_df <- all_samples_clr_df %>% 
  mutate(Residuals = resid(venom_lm))

# Create a second mi_df so that I don't have to add the following filtering step to all of the next ones
mi_df2 <- mi_df %>% 
  filter(In.Library == 'Yes') # Remove any proteins not in the library sent to Anthony
  


# Create a dataframe that has miRNAs and all the protein information in it that I can fuse later
miRNA_protein_pivot_df <- mi_df2 %>%
  dplyr::select(-Origin) %>% 
  distinct() %>% 
  filter(str_starts(Genes, 'Venom_')) %>%
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  grepl('VEGF', Genes) ~ 'VEGF',
                                  grepl('ohanin', Genes) ~ 'ohanin',
                                  grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('ADAM', Genes) ~ 'ADAM',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  grepl('CTL', Genes) ~ 'CTL',
                                  grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>% 
  dplyr::select(-contains('miRNA.Counts.'), -contains('RNA.VST'), -contains('RNA.Raw.Counts')) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('Intensity.'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>% 
  mutate(Sample.ID = gsub('Intensity.', '', Sample.ID)) %>% 
  distinct()

dup_rows1 <- miRNA_protein_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, Sample.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1)


# Create dataframe that has miRNAs in it that I can fuse to this one.
miRNA_mRNA_pivot_df <- mi_df2 %>% 
  dplyr::select(-Origin) %>% 
  distinct() %>% 
  filter(str_starts(Genes, 'Venom_')) %>%
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  grepl('VEGF', Genes) ~ 'VEGF',
                                  grepl('ohanin', Genes) ~ 'ohanin',
                                  grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('ADAM', Genes) ~ 'ADAM',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  grepl('CTL', Genes) ~ 'CTL',
                                  grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>% 
  dplyr::select(-contains('miRNA.Counts.'), -contains('Intensity'), -contains('RNA.Raw.Counts')) %>% 
  pivot_longer(
    cols = contains('RNA.VST'),
    names_to = 'Sample.ID',
    values_to = 'RNA.VST'
  ) %>% 
  mutate(Sample.ID = gsub('RNA.VST.', '', Sample.ID)) %>% 
  distinct()

# Identify duplicates
dup_rows2 <- miRNA_mRNA_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Start, miRNA.End, miRNA.Target.Start, Sample.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1)

# Create a variable to share common names between columns in the protein and mRNA data sets
shared_columns = intersect(names(miRNA_protein_pivot_df), names(miRNA_mRNA_pivot_df))

# Create a dataframe by fusing the protein and mRNA together
mRNA_protein_pivot_df <- left_join(
  miRNA_protein_pivot_df,
  miRNA_mRNA_pivot_df,
  by = shared_columns
) %>% 
  distinct()

# Identify duplicates
dup_rows3 <- mRNA_protein_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, Sample.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1)



# Create dataframe of only miRNAs and genes with their samples id
miRNA_pivot_df <- mi_df2 %>% 
  dplyr::select(-Origin) %>% 
  distinct() %>% 
  filter(str_starts(Genes, 'Venom_')) %>%
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  grepl('VEGF', Genes) ~ 'VEGF',
                                  grepl('ohanin', Genes) ~ 'ohanin',
                                  grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('ADAM', Genes) ~ 'ADAM',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  grepl('CTL', Genes) ~ 'CTL',
                                  grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>% 
  dplyr::select(-contains('miRNA.Counts.'), -contains('Intensity'), -contains('RNA.VST'), -contains('RNA.Raw.Counts')) %>% 
  pivot_longer(
    cols = contains('miRNA.RPM'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.RPM'
  ) %>%
  mutate(Sample.ID = gsub('miRNA.RPM.', '', Sample.ID)) %>% 
  distinct()

# Check for duplicated rows  
dup_rows4 <- miRNA_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, miRNA.Target.End, Venom.Family) %>%
  summarise(n = n()) %>%
  filter(n > 1)

# Get common columns between the two dataframes so they can be joined efficiently.
common_columns <- intersect(names(mRNA_protein_pivot_df), names(miRNA_pivot_df))

# Fuse dataframes back together
miRNA_mRNA_protein_pivot_df <- left_join(
  mRNA_protein_pivot_df,
  miRNA_pivot_df,
  by = common_columns
) %>% 
  dplyr::select(-contains('miRNA.RPM.')) %>% 
  distinct()

# Check for duplicated rows  
dup_rows5 <- miRNA_mRNA_protein_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, miRNA.Target.End, Venom.Family, Intensity, RNA.VST) %>%
  summarise(n = n()) %>%
  filter(n > 1)


# Create a set of shared columns
columns_shared <- intersect(names(miRNA_mRNA_protein_pivot_df), names(all_samples_residuals_df))

# Remove the following characters from the Sample.ID. column
characters_to_remove <- c(
  'LVG.2.',
  'LVG.9.',
  'LVG.4.',
  'RVG.5S.',
  'RVG.6S.',
  'RVG.7S.'
)
characters_to_remove2 <- c(
  '.viridis.Mid.M',
  '.viridis.North.M',
  '.viridis.South.M',
  '.viridis.North.F',
  '.concolor.Other.F',
  '.lutosus.Other.M'
)

# IMPORTANT: DON'T DELETE
# Fuse to the residuals dataframe
miRNA_with_residuals_df <- left_join(
  all_samples_residuals_df,
  miRNA_mRNA_protein_pivot_df,
  by = columns_shared
) %>% 
  distinct() %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove, collapse = "|"), '')) %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove2, collapse = "|"), ''))

# IMPORTANT: DON'T DELETE
# Create a wider version of the above
wider_miRNA_with_residuals_df <- miRNA_with_residuals_df %>% 
  pivot_wider(
    names_from = 'miRNA.Cluster',
    values_from = 'miRNA.RPM'
  ) %>% 
  dplyr::select(
    -contains('miRNA'), -contains('Blast'), -contains('Score'), -contains('Energy'), -contains('Scaled.'),
    -Positions, -Strand, -E.value, -RNA.VST, -Intensity, -In.Library, -Protein.Observed
  ) %>% 
  distinct()




#### Residuals Graphs (based on VST) ####
# 
# # Create a dataframe just for ohanin and filter at any missing data
# ohanin_df <- wider_miRNA_with_residuals_df %>% 
#   filter(Genes == 'Venom_ohanin') %>% 
#   dplyr::select(where(~any(!is.na(.))))
# 
# # Create plot for cvi-miR-145-5p
# # Create variable to control x and y of the ggplot
# miRNA_expression <- as.numeric(ohanin_df$Cluster_2576)
# exp_residuals <- abs(as.numeric(ohanin_df$Residuals))
# gene <- ohanin_df$Genes
# 
# # Create graph of values
# miR_145_5p_plot <- ohanin_df %>% 
#   ggplot(aes(x = miRNA_expression, y = exp_residuals)) +
#   geom_point(aes(color = gene), size = 2.5, alpha = 0.8) +
#   ggtitle('cvi-miR-145-5p miRNA Abundance compared to residuals') +
#   geom_smooth(
#     method = 'lm',
#     se = T,
#     color = 'black',
#     linetype = 'dashed',
#     formula = y ~ x - 1
#   ) +
#   stat_poly_eq(
#     aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
#     formula = y ~ x - 1
#   ) + 
#   xlab('miRNA Abundance (RPM)') +
#   ylab('| Residuals |') +
#   scale_color_manual(values = venom_colors) +
#   theme_classic() +
#   theme(legend.position = 'bottom',
#         legend.title = element_blank())
# miR_145_5p_plot
# 
# 


source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/Residuals_Correlation_Function.R')

# Create a list of venom genes for a function to iterate through
venom_genes <- wider_miRNA_with_residuals_df %>% distinct(Genes) 
venom_genes <- venom_genes$Genes
# Create path for the plots to go in
path <- 'Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/x_Individual_miRNA_Gene_Relationships'

# # Loop the function 
# for (gene in venom_genes) {
  # # Run the function in a loop
  # plots <- residuals_vs_miRNA_plot2(wider_miRNA_with_residuals_df, gene, path, '2024.09.14', filter_r_squared = F)
# }


# Add only the plots with R squared higher than 
# Create a second path
path2 <- 'Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/x_Individual_miRNA_Gene_Relationships/Good_R_squared'

# # Loop through the function
# for (gene in venom_genes) {
  # # Run the function in a loop
  # plots <- residuals_vs_miRNA_plot2(wider_miRNA_with_residuals_df, gene, path2, '2024.09.14', filter_r_squared = T)
# }



# # Define a helper function to call residuals_vs_miRNA_plot2 with different arguments
# run_plot <- function(gene, filter_r_squared, path) {
#   residuals_vs_miRNA_plot2(
#     data = wider_miRNA_with_residuals_df,
#     Gene = gene,
#     parent_directory = path,
#     Date = '2024.09.14',
#     filter_r_squared = filter_r_squared
#   )
# }

# # Run the function in parallel for the first case (no R-squared filtering)
# mclapply(venom_genes, run_plot, filter_r_squared = FALSE, path = path, mc.cores = num_cores)

# Define your batch size
batch_size <- 5  # Adjust this depending on your system's memory limits

# Split genes into batches
gene_batches <- split(venom_genes, ceiling(seq_along(venom_genes)/batch_size))

# Parallelize the execution of batches
mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot2(wider_miRNA_with_residuals_df, gene, path, '2024.09.14', filter_r_squared = FALSE)
    gc()
  }
}, mc.cores = num_cores)

# # Run the function in parallel for the second case (filtered by R-squared)
# mclapply(venom_genes, run_plot, filter_r_squared = TRUE, path = path2, mc.cores = num_cores)

mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot2(wider_miRNA_with_residuals_df, gene, path2, '2024.09.14', filter_r_squared = TRUE)
    gc()
  }
}, mc.cores = num_cores)





# #### Bubble Plot of miRNAs and Transcript-Protein Residuals ####
# 
# # Creata a data frame that only contains miRNA and Gene relationships for future left_join
# gene_mirna_relationships <- mi_df2 %>% 
#   dplyr::select(
#     Genes, miRNA.Cluster, Origin
#   ) %>% 
#   distinct()
# 
# # Import data frame with residuals and calculate variance for residuals
# residuals_variance_df <- miRNA_with_residuals_df %>% 
#   dplyr::select(
#     Genes, Sample.ID, Residuals
#   ) %>% 
#   group_by(Genes) %>% 
#   summarize(Residual.Variance = var(Residuals))
# 
# # Import data frame with miRNA expression levels and calculate variance
# mirna_variance_df <- miRNA_with_residuals_df %>% 
#   dplyr::select(
#     miRNA.Cluster, Sample.ID, miRNA.RPM
#   ) %>% 
#   group_by(miRNA.Cluster) %>% 
#   summarize(miRNA.RPM.Variance = var(miRNA.RPM))
# 
# # Fuse both variance data frames
# variance_residuals_mirna_df <- left_join(
#   residuals_variance_df,
#   gene_mirna_relationships,
#   by = c('Genes')
# ) %>% 
#   left_join(
#     mirna_variance_df,
#     by = c('miRNA.Cluster')
#   )
# 
# 
# # Create bubble plot that encodes residual variance with residual variance as bubble size and miRNA RPM variance as color
# residuals_mirna_expression_bubble_plot <- ggplot(variance_residuals_mirna_df, aes(x = Genes, y = miRNA.Cluster)) +
#   geom_point(aes(size = Residual.Variance, fill = log(miRNA.RPM.Variance + 1)), alpha = 0.75, shape = 21) +  # 'stroke' controls the width of the outline
#   scale_fill_viridis_c(option = 'viridis') +
#   scale_size_continuous(range = c(3, 13)) +
#   labs(
#     x = 'Genes',
#     y = 'miRNAs',
#     size = 'Residual Variance',
#     fill = 'miRNA Expression Variance',
#     title = 'Residuals and miRNA Expression Variance'
#   ) +
#   theme_classic2() +
#   theme(
#     plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
#     legend.key=element_blank(), 
#     axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
#     axis.text.y = element_text(colour = "black", size = 8), 
#     legend.text = element_text(size = 10, face ="bold", colour ="black"), 
#     legend.title = element_text(size = 12, face = "bold"), 
#     legend.position = "right"
#   )
# residuals_mirna_expression_bubble_plot
# ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Residuals_and_miRNA_Expression_2024.09.14.pdf", plot = residuals_mirna_expression_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)


#### Bubble Plot of R^2 values and Binding Score ####

# Create a data frame that contains the origin information
origin_df <- mi_df2 %>% 
  dplyr::select(
    miRNA.Cluster, miRNA.Start, miRNA.End, miRNA.Strandedness, miRNA.Target.Start, miRNA.Target.End, Total.Score, Total.Energy, Max.Score, Max.Energy, Origin
  )

# Create a shared set of names
shared_column_names_for_origin_fusion = intersect(names(origin_df), names(miRNA_with_residuals_df))


# Create a data frame that contains residauls and other important data
miRNA_with_residuals_and_origin_df <- left_join(
  miRNA_with_residuals_df,
  origin_df,
  by = shared_column_names_for_origin_fusion,
  relationship = 'many-to-many'
) %>% 
  dplyr::select(
    Sample.ID, Genes, Venom.Family, Residuals, miRNA.Cluster, Total.Score, Total.Energy, miRNA.RPM, Origin
  ) %>% 
  mutate(
    Venom.Color = case_when(
      grepl('SVMP', Genes) ~ SVMP_color,
      grepl('VEGF', Genes) ~ VEGF_color,
      grepl('ohanin', Genes) ~ ohanin_color,
      grepl('vQC', Genes) ~ vQC_color,
      grepl('SVSP', Genes) ~ SVSP_color,
      grepl('PLA2', Genes) ~ PLA2_color,
      grepl('CRISP', Genes) ~ CRISP_color,
      grepl('CTL', Genes) ~ CTL_color,
      grepl('EXO', Genes) ~ EXO_color,
      grepl('LAAO', Genes) ~ LAAO_color,
      grepl('myotoxin', Genes) ~ myotoxin_color,
      grepl('BPP', Genes) ~ BPP_color,
      TRUE ~ other_color
      )
  ) %>% 
  mutate(Genes = str_remove(Genes, '^Venom_'))

# Define a function to calculate the R squared value for each miRNA-Gene pair
calc_r_squared_and_slope <- function(data) {
  model <- lm(Residuals ~ miRNA.RPM, data = data)
  r_squared <- summary(model)$r.squared
  slope <- coef(model)['miRNA.RPM'] # Extract the slope (coefficient for miRNA.RPM)
  
  return(list(r_squared = r_squared, slope = slope))
}

# Calculate the R^2 values grouped by Genes and miRNA.Cluster
r_squared_df <- miRNA_with_residuals_and_origin_df %>% 
  group_by(Genes, miRNA.Cluster) %>% 
  summarize(result = list(calc_r_squared_and_slope(pick(everything()))), .groups = 'drop') %>%
  unnest_wider(result, names_sep = "_") %>%  # Unpack the list into separate columns
  rename(
    'R.Squared' = 'result_r_squared',
    'Slope' = 'result_slope'
  )

# Merge the R^2 values back into the original df
r_squared_df2 <- left_join(
  miRNA_with_residuals_and_origin_df,
  r_squared_df,
  by = c('Genes', 'miRNA.Cluster'),
) %>% 
  # dplyr::filter(
  #   miRNA.Cluster == c('cvi-miR-148a-3p', 'cvi-miR-737-5p'),
  #   Genes == 'Venom_CRISP3'
  # ) %>%  # This just checks that the R squared values are the same for the same miRNA-Gene relationship
  dplyr::select(
    -Sample.ID, -Residuals, -miRNA.RPM 
  ) %>% # This just removes unneccessary columns that could create repeats
  distinct() %>% 
  mutate(
    Origin.Color = case_when(
      grepl('three_prime_utr', Origin) ~ three_prime_color,
      grepl('five_prime_utr', Origin) ~ five_prime_color,
      grepl('CDS', Origin) ~ cds_color,
      TRUE ~ other_color
    )
  )


# View the updated data frame
print(head(r_squared_df2, n = 10))



# First lets try and see if binding score and binding energy correlate with the R^2 values at all
# Binding energy
be_vs_r_squared_plot <- ggplot(r_squared_df2, aes(x = abs(Total.Energy), y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'R Squared',
    title = 'Binding Energy vs R Squared'
  )
be_vs_r_squared_plot
# They do not correlate at all
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Binding_Data_vs_RSquared/miRNA_Binding_Energy_vs_RSquared_Regressions_2024.09.14.pdf", plot = be_vs_r_squared_plot, width = 10, height = 10, dpi = 900, create.dir = T)

# Create the same as above but include lables
be_vs_r_squared_plot2 <- ggplot(r_squared_df2, aes(x = abs(Total.Energy), y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  geom_text(
    aes(label = paste(Genes, miRNA.Cluster, sep = ": ")),  # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'R Squared',
    title = 'Binding Energy vs R Squared'
  )
be_vs_r_squared_plot2
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Binding_Data_vs_RSquared/Labeled/Labeled_miRNA_Binding_Energy_vs_RSquared_Regressions_2024.09.14.pdf", plot = be_vs_r_squared_plot2, width = 25, height = 25, dpi = 900, limitsize = F, create.dir = T)



# Binding Score
bs_vs_r_squared_plot <- ggplot(r_squared_df2, aes(x = Total.Score, y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'R Squared',
    title = 'Binding Score vs R Squared'
  )
bs_vs_r_squared_plot
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Binding_Data_vs_RSquared/miRNA_Binding_Score_vs_RSquared_Regressions_2024.09.14.pdf", plot = bs_vs_r_squared_plot, width = 10, height = 10, dpi = 900, create.dir = T)

# Same as above, but labeled to hell
bs_vs_r_squared_plot2 <- ggplot(r_squared_df2, aes(x = Total.Score, y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  geom_text(
    aes(label = paste(Genes, miRNA.Cluster, sep = ": ")),  # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'R Squared',
    title = 'Binding Score vs R Squared'
  )
bs_vs_r_squared_plot2
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Binding_Data_vs_RSquared/Labeled/Labeled_miRNA_Binding_Score_vs_RSquared_Regressions_2024.09.14.pdf", plot = bs_vs_r_squared_plot2, width = 25, height = 25, dpi = 900, create.dir = T)


# Rescale the R squared values so that negatively sloped miRNA gene relationships are negative
r_squared_df3 <- r_squared_df2 %>% 
  mutate(
    R.Squared = ifelse(Slope < 0, -R.Squared, R.Squared) # This adjust the sign of the R.Squared based on slope
  )


# Create bubble plot that encodes R^2 as the color and Binding Score as the size
r_squared_binding_score_bubble_plot <- ggplot(r_squared_df3, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = Total.Score, fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = Total.Score, fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0, 
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +  # Use custom colors
  scale_size_continuous(range = c(1, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'R Squared',
    title = 'R Squared and Binding Score'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/RSquared_and_Binding_Score_2024.09.13.pdf", plot = r_squared_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = T)


# Create bubble plot that encodes R^2 as the color and Binding Score as the size
r_squared_binding_energy_bubble_plot <- ggplot(r_squared_df3, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = abs(Total.Energy), fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = abs(Total.Energy), fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0, 
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +  # Use custom colors
  scale_size_continuous(range = c(1, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'R Squared',
    title = 'R Squared and Binding Energy'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +  
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_energy_bubble_plot
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/RSquared_and_Binding_Energy_2024.09.13.pdf", plot = r_squared_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = T)

#### Come back here bitch ####

# Set limits for grey and red color scale
min_r_squared = 0.75
max_total_energy = -7
min_total_score = 155

# Lets try the above again, but this time filtering by R^2 >= 0.5, and Binding Energy <= -7
filtered_r_squared_binding_energy_bubble_plot <- r_squared_df3 %>% 
  dplyr::filter(
    abs(R.Squared) >= min_r_squared,
    Total.Energy <= max_total_energy
  ) %>% 
  # dplyr::filter(
  #   R.Squared >= min_r_squared,
  #   Total.Energy <= max_total_energy
  # ) %>% 
  ggplot(aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = abs(Total.Energy), fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = abs(Total.Energy), fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_gradient(low = 'grey', high = 'red', limits = c(min_r_squared, 1)) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +  # Use custom colors
  scale_size_continuous(range = c(5, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'R Squared',
    title = 'R Squared and Binding Energy'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filtered_r_squared_binding_energy_bubble_plot
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Filtered_RSquared_and_Binding_Energy_2024.09.13.pdf", plot = filtered_r_squared_binding_energy_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)


# Let's do the same but filter out low binding energies instead
# Check the range of sizes in the filtered data
filtered_r_squared_df <- r_squared_df3 %>% 
  dplyr::filter(
    abs(R.Squared) >= min_r_squared,
    Total.Score >= min_total_score
  )

range_size <- range(filtered_r_squared_df$Total.Score, na.rm = TRUE)

filtered_r_squared_binding_score_bubble_plot <- ggplot(filtered_r_squared_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = Total.Score, fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +
  # geom_point(aes(size = Total.Score, fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +
  # scale_fill_viridis_c(option = 'magma') +
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_gradient(low = 'grey', high = 'red', limits = c(min_r_squared, 1)) +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +
  scale_size_continuous(
    range = c(1, 15),  # Adjust these values as needed
    breaks = seq(from = min(range_size), to = max(range_size), length.out = 5)  # Adjust breaks to include the size range
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'R Squared',
    title = 'R Squared and Binding Score'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filtered_r_squared_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Filtered_RSquared_and_Binding_Score_2024.09.13.pdf", plot = filtered_r_squared_binding_score_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)




# Lets do something else, lets use a "contact" heat map to show well Gene-miRNA pairs with strong evidence of an effect on protein expression
# Set binding score and energy
max_total_energy2 <- -20
min_total_score2 <- 160
min_r_squared2 = 0.5


# Create plot
R_squared_heatmap <- r_squared_df3 %>% 
  dplyr::filter(
    Total.Score >= min_total_score2,
    Total.Energy <= max_total_energy2,
    abs(R.Squared) >= min_r_squared2
  ) %>% # Filter by binding energy and score
  ggplot(aes(x = Genes, y = miRNA.Cluster)) +
  geom_tile(aes(fill = R.Squared)) +
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    fill = 'R Squared',
    title = 'R Squared Heatmap'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
R_squared_heatmap
ggsave("Figures/Expression_Plots/CLR_Transformed/Heatmaps/Filtered_RSquared_Heatmap_(BS=160_BE=-20)_2024.09.13.pdf", plot = R_squared_heatmap, width = 10, height = 15, dpi = 900, create.dir = T)


# Let's try maximizing the scores for this and redoing the figure
max_scores_df <- r_squared_df3 %>% 
  mutate(
    Total.Energy.Scaled = rescale(-Total.Energy), # Negate if lower energy is better
    Total.Score.Scaled = rescale(Total.Score)     # Normalize Total.Score between 0 and 1
  ) %>%
  # Combine the scaled scores into a ranking metric (here we use the sum)
  mutate(
    Total.Support = abs(R.Squared) + Total.Energy.Scaled + Total.Score.Scaled
  ) %>%
  # Sort by the combined score and select the top 30
  arrange(desc(Total.Support)) %>%
  slice_head(n = 30)
# View the top 10 genes
max_scores_df

# Remove unneccessary columns and save the data
max_scores_df2 <- max_scores_df %>% 
  dplyr::select(
    -contains('Color'), -Total.Energy.Scaled, -Total.Score.Scaled
  ) %>% 
  dplyr::rename(
    'Target.Type' = 'Origin'
  ) %>% 
  mutate(
    R.Squared = abs(R.Squared)
  )
write.csv(max_scores_df2, 'Tables/Expression_Tables/CLR/Top30_Best_miRNA-Gene_Pairs_2024.09.13.csv', row.names = F)

# Create plot
R_squared_heatmap_best_30 <- ggplot(max_scores_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_tile(aes(fill = R.Squared)) +
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    fill = 'R Squared',
    title = 'R Squared Heatmap'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
R_squared_heatmap_best_30
ggsave("Figures/Expression_Plots/CLR_Transformed/Heatmaps/Filtered_RSquared_Heatmap_Best_30_2024.09.13.pdf", plot = R_squared_heatmap_best_30, width = 10, height = 15, dpi = 900, create.dir = T)



# Let's try bubble plots with the same data
# Binding score
r_squared_binding_score_bubble_plot_top_30 <- ggplot(max_scores_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = Total.Score, fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +
  # geom_point(aes(size = Total.Score, fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +
  # scale_fill_viridis_c(option = 'magma') +
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +
  scale_size_continuous(
    range = c(1, 15)  # Adjust these values as needed
    ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'R Squared',
    title = 'R Squared and Binding Score'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_score_bubble_plot_top_30
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Top30_RSquared_and_Binding_Score_2024.09.13.pdf", plot = r_squared_binding_score_bubble_plot_top_30, width = 10, height = 15, dpi = 900, create.dir = T)


# Binding energy
r_squared_binding_energy_bubble_plot_top_30 <- ggplot(max_scores_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = abs(Total.Energy), fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +
  # geom_point(aes(size = abs(Total.Energy), fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +
  # scale_fill_viridis_c(option = 'magma') +
  scale_fill_gradient2(
    low = 'red', 
    mid = 'white', 
    high = 'blue', 
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +
  scale_size_continuous(
    range = c(1, 15)  # Adjust these values as needed
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'R Squared',
    title = 'R Squared and Binding Energy'
  ) +
  # theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_energy_bubble_plot_top_30
ggsave("Figures/Expression_Plots/CLR_Transformed/Bubble_Plot/Top30_RSquared_and_Binding_Energy_2024.09.13.pdf", plot = r_squared_binding_energy_bubble_plot_top_30, width = 10, height = 15, dpi = 900, create.dir = T)




#### mRNA vs Protein for all C. viridis ####

# Create new data frame for only the viridis samples
viridis_plot_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, '.viridis.'))


# Set variables for protein and RNA expression levels
viridis_mRNA_expression <- viridis_plot_df$CLR.Scaled.mRNA.Expression
viridis_protein_expression <- viridis_plot_df$CLR.Scaled.Protein.Expression
venom_family <- viridis_plot_df$Venom.Family


# Create figure for mRNA and Protein expression
viridis_only_plot <- viridis_plot_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
viridis_only_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/Viridis_Only_Individuals_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = viridis_only_plot, create.dir = T)



#### mRNA vs Protein for CV1081 (viridis) ####

# Create new data frame for only the viridis samples
CV1081_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, 'LVG.2.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV1081_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV1081_df$CLR.Scaled.Protein.Expression
venom_family <- CV1081_df$Venom.Family


# Create figure for mRNA and Protein expression
CV1081_plot <- CV1081_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV1081_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/CV1081_viridis_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = CV1081_plot, create.dir = T)


#### mRNA vs Protein for CV0857 (viridis) ####

# Create new data frame for only the viridis samples
CV0857_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, 'LVG.4.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV0857_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV0857_df$CLR.Scaled.Protein.Expression
venom_family <- CV0857_df$Venom.Family


# Create figure for mRNA and Protein expression
CV0857_plot <- CV0857_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV0857_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/CV0857_viridis_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = CV0857_plot, create.dir = T)



#### mRNA vs Protein for CV1086 (viridis) ####

# Create new data frame for only the viridis samples
CV1086_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, 'LVG.9.')) 


# Set variables for protein and RNA expression levels
mRNA_expression <- CV1086_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV1086_df$CLR.Scaled.Protein.Expression
venom_family <- CV1086_df$Venom.Family


# Create figure for mRNA and Protein expression
CV1086_plot <- CV1086_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV1086_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/CV1086_viridis_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = CV1086_plot, create.dir = T)



#### mRNA vs Protein for CV1087 (viridis) ####

# Create new data frame for only the viridis samples
CV1087_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, 'RVG.5S.')) 


# Set variables for protein and RNA expression levels
mRNA_expression <- CV1087_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV1087_df$CLR.Scaled.Protein.Expression
venom_family <- CV1087_df$Venom.Family


# Create figure for mRNA and Protein expression
CV1087_plot <- CV1087_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV1087_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/CV1087_viridis_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = CV1087_plot, create.dir = T)


#### mRNA vs Protein for CV0987 (lutosus) ####

# Create new data frame for only the viridis samples
CV0987_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, 'RVG.6S.')) 


# Set variables for protein and RNA expression levels
mRNA_expression <- CV0987_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV0987_df$CLR.Scaled.Protein.Expression
venom_family <- CV0987_df$Venom.Family


# Create figure for mRNA and Protein expression
CV0987_plot <- CV0987_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV0987_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/CV0987_lutosus_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = CV0987_plot, create.dir = T)



#### mRNA vs Protein for CV0985 (concolor) ####

# Create new data frame for only the viridis samples
CV0985_df <- all_samples_clr_df %>% 
  filter(str_detect(Sample.ID, 'RVG.7S.')) 


# Set variables for protein and RNA expression levels
mRNA_expression <- CV0985_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV0985_df$CLR.Scaled.Protein.Expression
venom_family <- CV0985_df$Venom.Family


# Create figure for mRNA and Protein expression
CV0985_plot <- CV0985_df %>% 
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) + 
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV0985_plot

# Save 
ggsave('Figures/Expression_Plots/CLR_Transformed/Linear_Regression_VST/CV0985_concolor_mRNAvsProtein_Expression_Plot_2024.09.14.pdf', plot = CV0985_plot, create.dir = T)


