# Last updated: 2025/02/13
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(stringr)

residuals_vs_miRNA_plot <- function(data, Gene, parent_directory, Date) {
  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)
  require(stringr)
  
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
  
  # Create color scheme for the venom genes
  venom_colors <- c(
    SVMP = SVMP_color,
    ADAM = ADAM_color,
    SVSP = SVSP_color,
    PLA2 = PLA2_color,
    miRNA = miRNA_color,
    VEGF = VEGF_color,
    ohanin = ohanin_color,
    myotoxin = myotoxin_color,
    vQC = vQC_color,
    CRISP = CRISP_color,
    CTL = CTL_color,
    EXO = EXO_color,
    LAAO = LAAO_color,
    BPP = BPP_color,
    others = other_color
  )
  
  
  # Filter the data frame to include only the specified Gene and columns starting with "cvi-" or "Cluster_"
  gene_df <- data %>%
    filter(Genes == Gene) %>%
    select(where(~any(!is.na(.)))) %>% 
    select(matches("^cvi-|^Cluster_"), residuals, Genes, Sample.ID, Venom.Family)
  
  # Extract the gene name
  Gene <- gene_df$Genes
  
  # Extract sample IDs
  sample_id <- gene_df$Sample.ID
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # Iterate over each column in the filtered data frame
  for (miRNA_cluster in colnames(gene_df)) {
    
    # Skip Genes and residuals columns
    if (miRNA_cluster %in% c("Genes", "residuals", 'Sample.ID', 'Venom.Family')) {
      next
    }
    
    # Extract the miRNA name from the column name
    miRNA_name <- str_replace(miRNA_cluster, "^cvi-", "")
    
    # Extract miRNA expression and residuals
    miRNA_expression <- as.numeric(gene_df[[miRNA_cluster]])
    exp_residuals <- as.numeric(gene_df$residuals)
    
    # Create the plot
    plot <- ggplot(data = gene_df, aes(x = miRNA_expression, y = exp_residuals, label = sample_id)) +
      geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
      geom_smooth(
        method = 'lm',
        se = TRUE,
        color = 'black',
        linetype = 'dashed',
        formula = y ~ x
        # formula = y ~ x - 1
      ) +
      stat_poly_eq(
        aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
        formula = y ~ x,
        # formula = y ~ x - 1,
        parse = TRUE
      ) + 
      geom_text(check_overlap = F, size = 3, hjust = -0.2) +
      labs(
        title = paste(miRNA_name, 'miRNA abundance vs residuals'),
        x = 'miRNA (rpm)',
        y = 'residuals'
      ) +
      scale_color_manual(values = venom_colors) +  # Apply color scheme
      theme_linedraw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, hjust = 0.5),
        legend.position = 'bottom',
        legend.title.position = 'top'
      )
    
    # Store the plot in the list
    plots[[miRNA_name]] <- plot
    
    # Date should be in the year.month.day format
    # Save the plot as a PDF file under the given parent_directory
    file_name <- paste(Gene, miRNA_name, "plot", Date, "pdf", sep = ".")
    ggsave(filename = file.path(parent_directory, file_name), plot = plot, device = "pdf", create.dir = T)
  }
  
  # Return the list of plots
  return(plots)
}


residuals_vs_miRNA_plot2 <- function(data, Gene, parent_directory, Date, filter_r_squared = FALSE) {
  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)
  require(stringr)
  
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
  
  # Create color scheme for the venom genes
  venom_colors <- c(
    SVMP = SVMP_color,
    ADAM = ADAM_color,
    SVSP = SVSP_color,
    PLA2 = PLA2_color,
    miRNA = miRNA_color,
    VEGF = VEGF_color,
    ohanin = ohanin_color,
    myotoxin = myotoxin_color,
    vQC = vQC_color,
    CRISP = CRISP_color,
    CTL = CTL_color,
    EXO = EXO_color,
    LAAO = LAAO_color,
    BPP = BPP_color,
    others = other_color
  )
  
  
  # Filter the data frame to include only the specified Gene and columns starting with "cvi-" or "Cluster_"
  gene_df <- data %>%
    filter(Genes == Gene) %>%
    select(where(~any(!is.na(.)))) %>% 
    select(matches("^cvi-|^Cluster_"), residuals, Genes, Sample.ID, Venom.Family)
  
  # Extract the gene name
  Gene <- gene_df$Genes
  
  # Extract sample IDs
  sample_id <- gene_df$Sample.ID
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # Iterate over each column in the filtered data frame
  for (miRNA_cluster in colnames(gene_df)) {
    
    # Skip Genes and residuals columns
    if (miRNA_cluster %in% c("Genes", "residuals", 'Sample.ID', 'Venom.Family')) {
      next
    }
    
    # Extract the miRNA name from the column name
    miRNA_name <- str_replace(miRNA_cluster, "^cvi-", "")
    
    # Extract miRNA expression and residuals
    miRNA_expression <- as.numeric(gene_df[[miRNA_cluster]])
    exp_residuals <- as.numeric(gene_df$residuals)
    
    
    # Perform linear regression so that only high correlation coefficients are used
    # Perform linear regression
    # lm_model <- lm(exp_residuals ~ miRNA_expression - 1)
    lm_model <- lm(exp_residuals ~ miRNA_expression)
    
    
    # Calculate R-squared value
    r_squared <- summary(lm_model)$r.squared
    
    # Skip plotting if filter_r_squared is TRUE and R-squared is less than 0.5
    if (filter_r_squared && r_squared < 0.5) {
      next
    }
    
    
    # Create the plot
    plot <- ggplot(data = gene_df, aes(x = miRNA_expression, y = exp_residuals, label = sample_id)) +
      geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
      geom_smooth(
        method = 'lm',
        se = TRUE,
        color = 'black',
        linetype = 'dashed',
        formula = y ~ x
        # formula = y ~ x - 1
      ) +
      stat_poly_eq(
        aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
        formula = y ~ x,
        # formula = y ~ x - 1,
        parse = TRUE
      ) + 
      geom_text(check_overlap = F, size =3, hjust = -0.2) +
      labs(
        title = paste(miRNA_name, 'miRNA abundance vs residuals'),
        x = 'miRNA (rpm)',
        y = 'residuals'
      ) +
      scale_color_manual(values = venom_colors) +  # Apply color scheme
      theme_linedraw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, hjust = 0.5),
        legend.position = 'bottom',
        legend.title.position = 'top'
      )
    
    # Store the plot in the list
    plots[[miRNA_name]] <- plot
    
    # Date should be in the year.month.day format
    # Save the plot as a PDF file under the given parent_directory
    file_name <- paste(Gene, miRNA_name, "plot", Date, "pdf", sep = ".")
    ggsave(filename = file.path(parent_directory, file_name), plot = plot, device = "pdf", create.dir = T)
  }
  
  # Return the list of plots
  return(plots)
}


residuals_vs_miRNA_plot3 <- function(data, gene, parent_directory, date, filter_r_squared = FALSE) {
  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)
  require(stringr)
  require(ggrepel)
  
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
  
  # Create color scheme for the venom genes
  venom_colors <- c(
    SVMP = SVMP_color,
    ADAM = ADAM_color,
    SVSP = SVSP_color,
    PLA2 = PLA2_color,
    miRNA = miRNA_color,
    VEGF = VEGF_color,
    ohanin = ohanin_color,
    myotoxin = myotoxin_color,
    vQC = vQC_color,
    CRISP = CRISP_color,
    CTL = CTL_color,
    EXO = EXO_color,
    LAAO = LAAO_color,
    BPP = BPP_color,
    others = other_color
  )
  
  
  # Filter the data frame to include only the specified data
  gene_df <- data %>%
    filter(genes == gene) %>%
    select(sample.id, genes, venom.family, residuals, miRNA.cluster, miRNA.rpm)
  
  # Extract the gene name
  gene <- gene_df$genes
  
  # Extract sample IDs
  sample_id <- gene_df$sample.id
  
  # Create a list of miRNA clusters that target the gene
  miRNA_clusters <- unique(gene_df$miRNA.cluster)
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # Iterate over each miRNA in the miRNA_clusters list
  for (miRNA_cluster in miRNA_clusters) {
    
    # Filter out everything in the gene_df other than the current miRNA of interest
    miRNA_df <- gene_df %>% 
      filter(miRNA.cluster == miRNA_cluster)
    
    # Extract miRNA expression and residuals
    miRNA_expression <- as.numeric(miRNA_df$miRNA.rpm)
    exp_residuals <- as.numeric(miRNA_df$residuals)
    
    # Perform linear regression
    lm_model <- lm(exp_residuals ~ miRNA_expression)
    
    # Calculate R-squared value
    r_squared <- summary(lm_model)$r.squared
    
    # Skip plotting if filter_r_squared is TRUE and R-squared is less than 0.5
    if (filter_r_squared && r_squared < 0.5) {
      next
    }
    
    # Create the plot
    plot <- ggplot(data = miRNA_df, aes(x = miRNA.rpm, y = residuals, label = sample.id, color = venom.family)) +
      geom_point(size = 2.5, alpha = 0.8) +
      # geom_point() +
      geom_smooth(
        method = 'lm',
        se = TRUE,
        linetype = 'dashed',
        formula = y ~ x
        # formula = y ~ x - 1
      ) +
      stat_poly_eq(
        aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
        formula = y ~ x,
        # formula = y ~ x - 1,
        parse = TRUE
      ) + 
      geom_text_repel(size = 3, hjust = -0.2) +
      labs(
        title = paste(miRNA_cluster, ' abundance vs residuals'),
        x = 'miRNA (rpm)',
        y = 'Protein Residuals',
        color = 'Venom Family'
      ) +
      scale_color_manual(values = venom_colors) +  # Apply color scheme
      theme_linedraw() +
      theme(
        plot.title = element_text(face = 'bold', size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, hjust = 0.5),
        legend.position = 'bottom',
        legend.title.position = 'top'
      )
    
    # Store the plot in the list
    plots[[miRNA_cluster]] <- plot
    
    # Date should be in the year.month.day format
    # Save the plot as a PDF file under the given parent_directory
    file_name <- paste(gene, miRNA_cluster, "plot", date, "pdf", sep = ".")
    ggsave(filename = file.path(parent_directory, file_name), plot = plot, device = "pdf", create.dir = T)
  }
  
  # Return the list of plots
  return(plots)
}



# Example usage:
# Assuming your wider_miRNA_with_residuals_df is named wider_miRNA_with_residuals_df
# Set the parent_directory where plots will be saved
# parent_directory <- "/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Figures/Expression_Plots/3UTR"
# 
# # Generate and save plots
# plots <- residuals_vs_miRNA_plot(wider_miRNA_with_residuals_df, 'Venom_ohanin', parent_directory, '2024.4.22')


# # Access plots for individual miRNAs
# plots$cvi_miR_1788_5p
# plots$cvi_another_miRNA