library(ggplot2)
library(dplyr)
library(ggpmisc)
library(stringr)

miRNA_vs_protein_plot <- function(data, Gene, parent_directory, Date) {
  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)
  require(stringr)
  
  # Filter the data frame to include only the specified Gene and columns starting with "cvi-" or "Cluster_"
  gene_df <- data %>%
    filter(Genes == Gene) %>%
    select(where(~any(!is.na(.)))) %>% 
    select(matches("^cvi-|^Cluster_"), Log.Intensity, Genes, Sample.ID)
  
  # Extract the gene name
  Gene <- gene_df$Genes
  
  # Extract sample IDs
  sample_id <- gene_df$Sample.ID
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # Iterate over each column in the filtered data frame
  for (miRNA_cluster in colnames(gene_df)) {
    
    # Skip Genes and Log.Intensity columns
    if (miRNA_cluster %in% c("Genes", "Log.Intensity", 'Sample.ID')) {
      next
    }
    
    # Extract the miRNA name from the column name
    miRNA_name <- str_replace(miRNA_cluster, "^cvi-", "")
    
    # Extract miRNA expression and residuals
    miRNA_expression <- as.numeric(gene_df[[miRNA_cluster]])
    Log.Intensity <- as.numeric(gene_df$Log.Intensity)
    
    # Create the plot
    plot <- ggplot(data = gene_df, aes(x = miRNA_expression, y = Log.Intensity, label = sample_id)) +
      geom_point(aes(color = Gene), size = 2.5, alpha = 0.8) +
      ggtitle(paste(miRNA_name, "miRNA abundance compared to mass spec Intensity (log scaled)")) +
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
      xlab('miRNA abundance (VST)') +
      ylab('Log.Intensity') +
      scale_color_brewer(palette = 'Set3') +
      theme_classic() +
      theme(legend.position = 'bottom',
            legend.title = element_blank())
    
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


miRNA_vs_protein_plot2 <- function(data, Gene, parent_directory, Date, filter_r_squared = FALSE) {
  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)
  require(stringr)
  
  # Filter the data frame to include only the specified Gene and columns starting with "cvi-" or "Cluster_"
  gene_df <- data %>%
    filter(Genes == Gene) %>%
    select(where(~any(!is.na(.)))) %>% 
    select(matches("^cvi-|^Cluster_"), Log.Intensity, Genes, Sample.ID)
  
  # Extract the gene name
  Gene <- gene_df$Genes
  
  # Extract sample IDs
  sample_id <- gene_df$Sample.ID
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # Iterate over each column in the filtered data frame
  for (miRNA_cluster in colnames(gene_df)) {
    
    # Skip Genes and Log.Intensity columns
    if (miRNA_cluster %in% c("Genes", "Log.Intensity", 'Sample.ID')) {
      next
    }
    
    # Extract the miRNA name from the column name
    miRNA_name <- str_replace(miRNA_cluster, "^cvi-", "")
    
    # Extract miRNA expression and residuals
    miRNA_expression <- as.numeric(gene_df[[miRNA_cluster]])
    Log.Intensity <- as.numeric(gene_df$Log.Intensity)
    
    
    # Perform linear regression so that only high correlation coefficients are used
    # Perform linear regression
    # lm_model <- lm(Log.Intensity ~ miRNA_expression - 1)
    lm_model <- lm(Log.Intensity ~ miRNA_expression)
    
    
    # Calculate R-squared value
    r_squared <- summary(lm_model)$r.squared
    
    # Skip plotting if filter_r_squared is TRUE and R-squared is less than 0.5
    if (filter_r_squared && r_squared < 0.5) {
      next
    }
    
    
    # Create the plot
    plot <- ggplot(data = gene_df, aes(x = miRNA_expression, y = Log.Intensity, label = sample_id)) +
      geom_point(aes(color = Gene), size = 2.5, alpha = 0.8) +
      ggtitle(paste(miRNA_name, "miRNA abundance compared to mass spec intensity (log scaled)")) +
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
      xlab('miRNA abundance (VST)') +
      ylab('Log.Intensity') +
      scale_color_brewer(palette = 'Set3') +
      theme_classic() +
      theme(legend.position = 'bottom',
            legend.title = element_blank())
    
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


# Example usage:
# Assuming your wider_miRNA_with_residuals_df is named wider_miRNA_with_residuals_df
# Set the parent_directory where plots will be saved
# parent_directory <- "/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Figures/Expression_Plots/3UTR"
# 
# # Generate and save plots
# plots <- miRNA_vs_protein_plot(wider_miRNA_with_residuals_df, 'Venom_ohanin', parent_directory, '2024.4.22')


# # Access plots for individual miRNAs
# plots$cvi_miR_1788_5p
# plots$cvi_another_miRNA