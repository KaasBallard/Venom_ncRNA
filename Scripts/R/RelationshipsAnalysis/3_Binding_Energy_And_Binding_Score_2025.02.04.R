# Last Edited: 2025/02/04

# Set up and Read Data ----

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
library(arrow)

# Create variable for the fused data set.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read in the data
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    !is.na(feature.type)
  ) %>% # This should get rid of all of the weirdly annotated genes
  distinct()
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
EXO_color <- '#005824'
LAAO_color <- '#B35806'
BPP_color <- '#1B9E77'
other_color <- '#666666'
three_prime_color <- 'black'
five_prime_color <- '#0072b2'
cds_color <- '#d55e00'

# Create an order for the genes
venom_gene_order <- c(
  'BPP',
  'myotoxin',
  'ohanin',
  'PLA2B1', 'PLA2K', 'PLA2C1', 'PLA2A1',
  'SVSP1', 'SVSP2', 'SVSP3', 'SVSP4', 'SVSP10', 'SVSP5', 'SVSP6', 'SVSP11', 'SVSP7', 'SVSP8', 'SVSP9',
  'SVMP1', 'SVMP2', 'SVMP3', 'SVMP4', 'SVMP5', 'SVMP6', 'SVMP7', 'SVMP8', 'SVMP9', 'SVMP10', 'SVMP11', 'SVMP12',
  'CRISP1', 'CRISP2', 'CRISP3', 'CRISP4',
  'CTL1', 'CTL2', 'CTL3', 'CTL4', 'CTL5', 'CTL6',
  'EXO1', 'EXO2', 'EXO3',
  'LAAO1', 'LAAO2', 'LAAO3',
  'VEGF1', 'VEGF2',
  'vQC1', 'vQC2'
) 


# Binding Energy vs Binding Score ----

# Create a data frame that just contains binding scores and energies, plus binding type 
be_bs_df <- mi_df %>% 
  filter(
    str_detect(genes, 'Venom_'),
    !str_detect(genes, 'ADAM')
  ) %>%
  mutate(color = case_when(
    grepl('SVMP', genes) ~ SVMP_color,
    grepl('ADAM', genes) ~ ADAM_color,
    grepl('VEGF', genes) ~ VEGF_color,
    grepl('ohanin', genes) ~ ohanin_color,
    grepl('vQC', genes) ~ vQC_color,
    grepl('SVSP', genes) ~ SVSP_color,
    grepl('PLA2', genes) ~ PLA2_color,
    grepl('CRISP', genes) ~ CRISP_color,
    grepl('CTL', genes) ~ CTL_color,
    grepl('EXO', genes) ~ EXO_color,
    grepl('LAAO', genes) ~ LAAO_color,
    grepl('myotoxin', genes) ~ myotoxin_color,
    grepl('BPP', genes) ~ BPP_color,
    TRUE ~ other_color
  )) %>% 
  dplyr::select(miRNA.cluster, genes, total.score, total.energy, feature.type, venom.family, color) %>% 
  distinct() %>% 
  mutate(
    genes = str_remove(genes,'^Venom_')
  )

# Enforce order for genes
be_bs_df$genes <- factor(be_bs_df$genes, levels = venom_gene_order)

# Create graph of total.energy as an explanation for total.score. Note that binding energy uses the absolute value to be more intuitive
be_vs_bs_plot <- ggplot(be_bs_df, aes(x = abs(total.energy), y = total.score)) +
  geom_point(aes(color = feature.type)) +
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
  scale_y_continuous(breaks = pretty(be_bs_df$total.score, n = 10)) +  # Increase y-axis labels
  theme_classic() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 12)
  ) +
  labs(
    x = 'Binding Energy',
    y = 'Binding Score',
    title = 'Binding Energy vs Binding Score'
  )
be_vs_bs_plot
ggsave("Figures/Binding_Score_Plots/Regressions/miRNA_Binding_Energy_vs_Binding_Score_Plot_2025.02.04.pdf", plot = be_vs_bs_plot, width = 10, height = 10, dpi = 900, create.dir = T)
# Wow, they are not super correlated and there are some clear distictions between types
# CDS are over represented in low quality scores
# There is also a lot of stuff on the low end of the quality spectrum and some high end stuff and not much in between

# Let's try the same thing, but remove some the CDS to see if the five prime utr is over repressented
be_vs_bs_plot_no_cds <- be_bs_df %>% 
  filter(!str_detect(feature.type, 'CDS')) %>% 
  ggplot(aes(x = abs(total.energy), y = total.score)) +
  geom_point(aes(color = feature.type)) +
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
  scale_y_continuous(breaks = pretty(be_bs_df$total.score, n = 10)) +  # Increase y-axis labels
  theme_classic() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 12)
  ) +
  labs(
    x = 'Binding Energy',
    y = 'Binding Score',
    title = 'Binding Energy vs Binding Score'
  )
be_vs_bs_plot_no_cds
# So it seems they are over represented



# Let's try adding venom family in for color instead
venom_family_be_bs_plot <- ggplot(be_bs_df, aes(x = abs(total.energy), y = total.score)) +
  geom_point(aes(color = venom.family, shape = feature.type)) +
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
  scale_color_manual(values = setNames(be_bs_df$color, be_bs_df$venom.family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(be_bs_df$total.score, n = 20)) +  # Increase y-axis labels
  theme_classic() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 12)
  ) +
  labs(
    x = 'Binding Energy',
    y = 'Binding Score',
    title = 'Binding Energy vs Binding Score'
  )
venom_family_be_bs_plot
# I don't see any pattern
ggsave("Figures/Binding_Score_Plots/Regressions/miRNA_Binding_Energy_vs_Binding_Score_Plot_Colored_by_Venom_Family_2025.02.04.pdf", plot = venom_family_be_bs_plot, width = 10, height = 10, dpi = 900, create.dir = T)


# Let's try splitting the regression lines between 250
split_be_bs_plot <- ggplot(be_bs_df, aes(x = abs(total.energy), y = total.score)) +
  geom_point(aes(color = venom.family, shape = feature.type)) +
  geom_smooth(
    aes(group = ifelse(total.score <= 250, 'Below 250', 'Above 250')), # Split the data based on Binding Score
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~"),
        group = ifelse(total.score <= 250, "Below 250", "Above 250")),  # Split the labels based on total.score
    formula = y ~ x
  ) +
  scale_color_manual(values = setNames(be_bs_df$color, be_bs_df$venom.family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(be_bs_df$total.score, n = 20)) +  # Increase y-axis labels
  theme_classic() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 12)
  ) +
  labs(
    x = 'Binding Energy',
    y = 'Binding Score',
    title = 'Binding Energy vs Binding Score'
  )
split_be_bs_plot
ggsave("Figures/Binding_Score_Plots/Regressions/miRNA_Binding_Energy_vs_Binding_Score_Plot_Split_Regressions_2025.02.04.pdf", plot = split_be_bs_plot, width = 10, height = 10, dpi = 900, create.dir = T)



# miRNA vs Gene Bubble Plot ----

# Calculate the ranges for each venom family
family_ranges <- be_bs_df %>%
  group_by(venom.family) %>%
  summarize(
    xmin = min(as.integer(genes)) - 0.3, # Convert factor levels to integers
    xmax = max(as.integer(genes)) + 0.3 # Add padding
  ) 
  
# Create a color vector
colors <- setNames(be_bs_df$color, be_bs_df$venom.family)

# Create a bubble plot with outline color based on miRNA target type
be_bs_bubble_plot <- ggplot(be_bs_df, aes(y = miRNA.cluster, x = factor(genes))) +
  geom_rect(
    data = family_ranges, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, color = venom.family),
    fill = 'black', alpha = 0.1, inherit.aes = FALSE
  ) + # Adjust transparency with alpha
  scale_color_manual(values = colors, name = 'Venom Family') +  # Use manual colors for venom.family
  geom_point(aes(size = total.score, fill = abs(total.energy)), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = total.score, fill = abs(total.energy), color = feature.type), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_viridis_c(option = 'magma') +
  scale_size_continuous(range = c(3, 13)) +
  labs(
    x = 'genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Binding Energy',
    title = 'Binding Energy and Binding Score'
  ) +
  theme_classic2() +
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
    feature.type ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      feature.type = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
be_bs_bubble_plot



# Try out different colors
magma <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'magma')
magma
inferno <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'inferno')
inferno
plasma <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'plasma')
plasma
viridis <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'viridis')
viridis
cividis <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'cividis')
cividis
rocket <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'rocket')
rocket
mako <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'mako')
mako
turbo <- be_bs_bubble_plot + scale_fill_viridis_c(option = 'turbo')
turbo
ggsave("Figures/Binding_Score_Plots/Bubble_Plots/miRNA_Binding_Energy_And_Binding_Score_Plot_2025.02.04.pdf", plot = be_bs_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)
# This is kind of to much, so lets try a bubble plot for each type of miRNA binding site


# Remove everything other than three prime utr targeting miRNAs
three_utr_be_bs_df <- be_bs_df %>% 
  filter(feature.type == 'three_prime_utr') %>% 
  distinct()

# Create a bubble plot for three prime targeting miRNAs
three_prime_be_bs_bubble_plot <- ggplot(three_utr_be_bs_df, aes(y = miRNA.cluster, x = genes)) +
  geom_point(aes(size = total.score, fill = abs(total.energy)), alpha = 0.75, shape = 21) +  # 'stroke' controls the width of the outline
  scale_fill_viridis_c(option = 'magma') +
  scale_size_continuous(range = c(3, 13)) +
  labs(
    x = 'genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Binding Energy',
    shape = 'feature.type',
    title = 'Binding Energy and Binding Score'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  )
three_prime_be_bs_bubble_plot
ggsave("Figures/Binding_Score_Plots/Bubble_Plots/three_prime_binding_miRNA_Binding_Energy_And_Binding_Score_Plot_2025.02.04.pdf", plot = three_prime_be_bs_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)



# Remove everything other than five prime utr targeting miRNAs
five_utr_be_bs_df <- be_bs_df %>% 
  filter(feature.type == 'five_prime_utr') %>% 
  distinct()

# Create a bubble plot for five prime targeting miRNAs
five_prime_be_bs_bubble_plot <- ggplot(five_utr_be_bs_df, aes(y = miRNA.cluster, x = genes)) +
  geom_point(aes(size = total.score, fill = abs(total.energy)), alpha = 0.75, shape = 21) +  # 'stroke' controls the width of the outline
  scale_fill_viridis_c(option = 'magma') +
  scale_size_continuous(range = c(3, 13)) +
  labs(
    x = 'genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Binding Energy',
    shape = 'feature.type',
    title = 'Binding Energy and Binding Score'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  )
five_prime_be_bs_bubble_plot
ggsave("Figures/Binding_Score_Plots/Bubble_Plots/five_prime_binding_miRNA_Binding_Energy_And_Binding_Score_Plot_2025.02.04.pdf", plot = five_prime_be_bs_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)



# Remove everything other than five prime utr targeting miRNAs
cds_be_bs_df <- be_bs_df %>% 
  filter(feature.type == 'CDS') %>% 
  distinct()

# Create a bubble plot for five prime targeting miRNAs
cds_be_bs_bubble_plot <- ggplot(cds_be_bs_df, aes(y = miRNA.cluster, x = genes)) +
  geom_point(aes(size = total.score, fill = abs(total.energy)), alpha = 0.75, shape = 21) +  # 'stroke' controls the width of the outline
  scale_fill_viridis_c(option = 'magma') +
  scale_size_continuous(range = c(3, 13)) +
  labs(
    x = 'genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Binding Energy',
    shape = 'feature.type',
    title = 'Binding Energy and Binding Score'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1), 
    axis.text.y = element_text(colour = "black", size = 8), 
    legend.text = element_text(size = 10, face ="bold", colour ="black"), 
    legend.title = element_text(size = 12, face = "bold"), 
    legend.position = "right"
  )
cds_be_bs_bubble_plot
ggsave("Figures/Binding_Score_Plots/Bubble_Plots/CDS_binding_miRNA_Binding_Energy_And_Binding_Score_Plot_2025.02.04.pdf", plot = cds_be_bs_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)
# Not sure what to go with for this stuff honestly




#### Filtered Binding Score VS Binding Energy Graph ####

# Filter the binding score and energy data frame by those numbers
# Set binding score and energy thresholds
binding_score <- 155
binding_energy <- -7

# Filter
filtered_be_bs_df <- be_bs_df %>% 
  filter(
    total.score >= binding_score,
    total.energy <= binding_energy
  )

# Update venom_gene_order to include only genes present in filtered data
updated_venom_gene_order <- venom_gene_order[venom_gene_order %in% filtered_be_bs_df$genes]

# Reapply the updated gene order to the filtered data
filtered_be_bs_df <- filtered_be_bs_df %>%
  mutate(genes = factor(genes, levels = updated_venom_gene_order))  # Use updated order

# Now recalculate the family ranges
family_ranges <- filtered_be_bs_df %>%
  group_by(venom.family) %>%
  summarize(
    xmin = min(as.integer(genes)) - 0.3,  # Convert factor levels to integers
    xmax = max(as.integer(genes)) + 0.3   # Add padding
  )

# Create a color vector
colors <- setNames(filtered_be_bs_df$color, filtered_be_bs_df$venom.family)

# Create a bubble plot with outline color based on miRNA target type
filtered_be_bs_bubble_plot <- ggplot(filtered_be_bs_df, aes(y = miRNA.cluster, x = factor(genes))) +
  geom_rect(
    data = family_ranges, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, color = venom.family),
    fill = 'black', alpha = 0.1, inherit.aes = FALSE
  ) + # Adjust transparency with alpha
  scale_color_manual(values = colors, name = 'Venom Family') +  # Use manual colors for venom.family
  geom_point(aes(size = total.score, fill = abs(total.energy)), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = total.score, fill = abs(total.energy), color = feature.type), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_viridis_c(option = 'magma') +
  scale_size_continuous(range = c(3, 13)) +
  labs(
    x = 'genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Binding Energy',
    title = 'b. Binding Energy and Binding Score'
  ) +
  theme_classic2() +
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
    feature.type ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      feature.type = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filtered_be_bs_bubble_plot
ggsave("Figures/Binding_Score_Plots/Bubble_Plots/Filtered_miRNA_Binding_Energy_And_Binding_Score_Plot_2025.02.04.pdf", plot = filtered_be_bs_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)
figure_e <- filtered_be_bs_bubble_plot
ggsave("Figures/1_Main_Figures/Figure_e/Figure_e_Bubble_Plot_2025.02.04.pdf", plot = figure_e, width = 10, height = 15, dpi = 900, create.dir = T)
