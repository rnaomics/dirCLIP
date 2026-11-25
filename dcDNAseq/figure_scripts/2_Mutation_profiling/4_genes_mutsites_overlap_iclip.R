#!/usr/bin/env Rscript
# Genes with Mutation Sites Overlapping iCLIP Peaks
# Manuscript-specific figure-generation script.
#
# Purpose:
# Compare the number of genes whose mutation sites overlap
# or do not overlap iCLIP peaks under different controls.

############################################################
# Number of genes with mutation sites overlapping iCLIP sites
############################################################

##############################
# Set working environment
##############################
# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")

##############################
# Set library path
##############################
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

##############################
# Load required packages
##############################
library(ggplot2)
library(scales)
library(dplyr)
library(grid)      # Required for unit()

############################################################
# Function: Create stacked bar plot
############################################################
create_stacked_bar_plot <- function(plot_df,
                                    fill_colors,
                                    plot_title,
                                    y_limits,
                                    y_tick_interval,
                                    x_label,
                                    y_label,
                                    legend_title) {
  
  text_size <- 16
  
  ##########################################################
  # Generate stacked bar plot
  ##########################################################
  figure <- ggplot(
    plot_df,
    aes(x = Group, y = y, fill = sub_group)
  ) +
    geom_bar(
      stat = "identity",
      width = 0.7,
      color = "black",
      linewidth = 0.6
    ) +
    labs(
      x = x_label,
      y = y_label,
      fill = legend_title,
      title = plot_title
    ) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal(base_size = text_size) +
    theme(
      # Remove grid lines
      panel.grid = element_blank(),
      
      # Draw only axis lines
      panel.border = element_blank(),
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      legend.position = "right",
      
      axis.title = element_text(size = text_size),
      axis.text = element_text(
        size = text_size,
        color = "black"
      ),
      
      legend.title = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      
      # Draw axis ticks
      axis.ticks = element_line(
        color = "black",
        size = 0.8
      ),
      
      # Tick length
      axis.ticks.length = unit(0.25, "cm"),
      
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black")
    )
  
  ##########################################################
  # Define y-axis limits and tick marks
  ##########################################################
  if (length(y_limits) > 0) {
    
    y_min <- y_limits[1]
    y_max <- y_limits[2]
    
  } else {
    
    group_totals <- plot_df %>%
      group_by(Group) %>%
      summarise(
        total_y = sum(y),
        .groups = "drop"
      )
    
    y_min <- 0
    
    y_max <- ceiling(
      max(group_totals$total_y) /
        y_tick_interval
    ) * y_tick_interval
  }
  
  y_ticks <- seq(
    y_min,
    y_max,
    by = y_tick_interval
  )
  
  figure <- figure +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(mult = c(0, 0.02))
    )
  
  return(figure)
}

############################################################
# Function:
# Plot genes containing mutation sites that overlap
# iCLIP sites
############################################################
plot_gene_overlap_with_iclip <- function(
    mutation_file_1,
    sample_label_1,
    mutation_file_2,
    sample_label_2,
    fill_colors,
    y_limits,
    y_tick_interval,
    peak_caller,
    output_file,
    figure_size) {
  
  ##########################################################
  # Load gene-level mutation overlap summary tables
  ##########################################################
  gene_overlap_df_1 <- read.delim(mutation_file_1)
  gene_overlap_df_2 <- read.delim(mutation_file_2)
  
  ##########################################################
  # Count genes with at least one mutation site that
  # overlaps an iCLIP site
  ##########################################################
  n_overlap_genes_1 <- sum(
    gene_overlap_df_1[
      ,
      paste0(
        "Sites_overlap_with_",
        peak_caller
      )
    ] > 0
  )
  
  n_overlap_genes_2 <- sum(
    gene_overlap_df_2[
      ,
      paste0(
        "Sites_overlap_with_",
        peak_caller
      )
    ] > 0
  )
  
  ##########################################################
  # Create counts for stacked bars
  ##########################################################
  gene_counts <- c(
    n_overlap_genes_1,
    sum(gene_overlap_df_1$N_sites > 0) -
      n_overlap_genes_1,
    
    n_overlap_genes_2,
    sum(gene_overlap_df_2$N_sites > 0) -
      n_overlap_genes_2
  )
  
  ##########################################################
  # Define sample groups
  ##########################################################
  group_labels <- c(
    sample_label_1,
    sample_label_1,
    sample_label_2,
    sample_label_2
  )
  
  ##########################################################
  # Define overlap categories
  ##########################################################
  overlap_status <- c(
    "Yes",
    "No",
    "Yes",
    "No"
  )
  
  ##########################################################
  # Build plotting data frame
  ##########################################################
  plot_df <- data.frame(
    y = gene_counts,
    Group = group_labels,
    sub_group = overlap_status
  )
  
  ##########################################################
  # Generate figure
  ##########################################################
  figure <- create_stacked_bar_plot(
    plot_df = plot_df,
    fill_colors = fill_colors,
    plot_title = "",
    y_limits = y_limits,
    y_tick_interval = y_tick_interval,
    x_label = "Background sample",
    y_label = "Number of genes",
    legend_title = "Overlap with iCLIP"
  )
  
  ##########################################################
  # Save figure
  ##########################################################
  ggsave(
    filename = output_file,
    plot = figure,
    width = figure_size[1],
    height = figure_size[2],
    limitsize = FALSE,
    device = cairo_pdf
  )
}

############################################################
# Input files
############################################################

# -UV sample
mutation_files_1 <- c(
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/cov_thresh_50/gencode/GW_mut_thresh_0.03.tsv"
)

sample_labels_1 <- c("-UV")

# +UV sample
mutation_files_2 <- c(
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_UV/2_Inc_del/cov_thresh_50/gencode/GW_mut_thresh_0.03.tsv"
)

sample_labels_2 <- c("+UV")

############################################################
# Output files
############################################################
output_files <- c(
  "Result/figure/1_Mutation_profiling/SRSF3_Genes_with_Mutsites_overlap_with_iCLIPsites_bar.pdf"
)

############################################################
# Analysis parameters
############################################################
peak_caller <- "iCount"

fill_colors <- c(
  "Yes" = "gray40",
  "No"  = "gray70"
)

############################################################
# Generate plots
############################################################
for (file_idx in seq_along(mutation_files_1)) {
  
  mutation_file_1 <- mutation_files_1[file_idx]
  mutation_file_2 <- mutation_files_2[file_idx]
  
  sample_label_1 <- sample_labels_1[file_idx]
  sample_label_2 <- sample_labels_2[file_idx]
  
  output_file <- output_files[file_idx]
  
  print(mutation_file_1)
  print(mutation_file_2)
  
  ##########################################################
  # Plot genes with mutation sites overlapping iCLIP sites
  ##########################################################
  y_tick_interval <- 500
  y_limits <- c()
  
  plot_gene_overlap_with_iclip(
    mutation_file_1 = mutation_file_1,
    sample_label_1 = sample_label_1,
    mutation_file_2 = mutation_file_2,
    sample_label_2 = sample_label_2,
    fill_colors = fill_colors,
    y_limits = y_limits,
    y_tick_interval = y_tick_interval,
    peak_caller = peak_caller,
    output_file = output_file,
    figure_size = c(6, 7)
  )
}

############################################################
# End of script
############################################################