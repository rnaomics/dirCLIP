################################################################################
# genes_with_mutsite_iclip_overlap.R
# Author: Mohamed B
# Description:
# For each peak caller, mutation type, and coverage threshold, calculates
# the number of genes containing mutation sites overlapping iCLIP peaks.
# Generates grouped bar plots showing overlap vs non-overlap distributions.
################################################################################

### 🧬 Required Libraries ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

library(Rsubread)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(GGally)
library(ggcorrplot)
library(scales)
library(RColorBrewer)
library(S4Vectors)
library(GenomicRanges)
library(writexl)

################################################################################
# 📊 Helper Function: Styled Bar Plot Generator
################################################################################

make_bar_plot <- function(df, Colors, Title, y_lim, y_tick_interval,
                          x_lab, y_lab, legend_title) {
  text_size <- 16
  fig <- ggplot(df, aes(x = Group, y = y, fill = sub_group)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    labs(x = x_lab, y = y_lab, fill = legend_title, title = Title) +
    scale_fill_manual(values = Colors) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(size = text_size),
      axis.text  = element_text(size = text_size, color = "black"),
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size),
      legend.position = "right",
      axis.ticks        = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    )
  
  ## Compute Y axis limits
  if (length(y_lim) == 0) {
    summary_df <- df %>% group_by(Group) %>% summarise(total_y = sum(y), .groups = "drop")
    y_min <- 0
    y_max <- ceiling(max(summary_df$total_y) / y_tick_interval) * y_tick_interval
  } else {
    y_min <- y_lim[1]
    y_max <- y_lim[2]
  }
  
  y_ticks <- seq(y_min, y_max, by = y_tick_interval)
  fig <- fig + scale_y_continuous(limits = c(y_min, y_max), breaks = y_ticks,
                                  expand = expansion(mult = c(0, 0.02)))
  return(fig)
}

###### 🧠 Main Plotting Function: Genes Overlapping iCLIP Sites ##############################################
plot_genes_with_mutsites_overlap_with_iclip_joined_plots <- function(Mut_path1, Mut_path2, filt_by,
                                                                     Coverage_thresholds, Diff_thresholds,
                                                                     Colors, y_lim, y_tick_interval,
                                                                     Peak_caller, fig_name, fig_sz, sub_plots) {
  fig <- ggdraw()
  sub_fig_sz <- 1 / sub_plots
  grid_pos <- expand.grid(x = 0:(sub_plots[1]-1), y = 0:(sub_plots[2]-1))
  
  for (k in seq_along(Coverage_thresholds)) {
    for (m in seq_along(Diff_thresholds)) {
      Cov_thresh <- Coverage_thresholds[k]
      Diff_thresh <- Diff_thresholds[m]
      
      # 📥 Load overlap tables
      file1 <- paste0(Mut_path1, "cov_thresh_", Cov_thresh, "/", filt_by, "/GW_mut_thresh_", Diff_thresh, ".tsv")
      file2 <- paste0(Mut_path2, "cov_thresh_", Cov_thresh, "/", filt_by, "/GW_mut_thresh_", Diff_thresh, ".tsv")
      Gene_olaps1 <- read.delim(file1)
      Gene_olaps2 <- read.delim(file2)
      
      # 🧮 Count genes with overlap
      N_Genes_1 <- colSums(Gene_olaps1[, -(1:5)] > 0)
      N_Genes_2 <- colSums(Gene_olaps2[, -(1:5)] > 0)
      
      key <- paste0("Sites_overlap_with_", Peak_caller)
      y_vals <- c(N_Genes_1[key], sum(Gene_olaps1$N_sites > 0) - N_Genes_1[key],
                  N_Genes_2[key], sum(Gene_olaps2$N_sites > 0) - N_Genes_2[key])
      
      df <- data.frame(
        y = y_vals,
        Group = c(Samples1, Samples1, Samples2, Samples2),
        sub_group = rep(c("Yes", "No"), 2)
      )
      
      Title <- paste0("Cov_thresh = ", Cov_thresh, ", Diff_thresh = ", Diff_thresh)
      bar <- make_bar_plot(df, Colors, Title, y_lim, y_tick_interval,
                           x_lab = "Background sample",
                           y_lab = "Number of genes",
                           legend_title = "Overlap with iCLIP")
      
      fig <- fig + draw_plot(bar,
                             x = grid_pos$x[k] * sub_fig_sz[1],
                             y = grid_pos$y[m] * sub_fig_sz[2],
                             width = sub_fig_sz[1],
                             height = sub_fig_sz[2])
    }
  }
  
  ## 💾 Save final figure
  ggsave(fig_name,
         plot = fig,
         width  = fig_sz[1] * sub_plots[1],
         height = fig_sz[2] * sub_plots[2],
         limitsize = FALSE,
         device = cairo_pdf)
}

##############################################################################################
## Input configuration
paths1         <- "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/"
ref_subfolder1 <- "ref_TotalRNA_cDNA_no_UV"
Samples1       <- "-UV"

paths2         <- "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/"
ref_subfolder2 <- "ref_TotalRNA_cDNA_UV"
Samples2       <- "+UV"

subfolders1 <- "q_8_Q_8"
subfolders2 <- c("2_Inc_del", "1_No_del")
Peak_callers <- c("CLIPper", "iCount", "pureclip", "piranha", "paraclu")

Coverage_thresholds <- c(50)
Diff_thresholds <- c(0.03)
sub_plots <- c(1, 1)
filt_by <- "FDR"

Colors <- c("Yes" = "gray40", "No" = "gray70")

## Loop over peak callers and conditions
for (Peak_caller in Peak_callers) {
  for (subfolder1 in subfolders1) {
    for (subfolder2 in subfolders2) {
      
      Mut_path1 <- paste0(paths1, subfolder1, "/", subfolder2, "/", ref_subfolder1, "/")
      Mut_path2 <- paste0(paths2, subfolder1, "/", subfolder2, "/", ref_subfolder2, "/")
      Res_path  <- paste0(Mut_path1, filt_by, "/")
      dir.create(Res_path, showWarnings = FALSE, recursive = TRUE)
      
      ## Plot bar chart
      fig_file <- paste0(Res_path, "Genes_with_Mutsites_overlap_with_", Peak_caller, "_iCLIPsites_bar.pdf")
      plot_genes_with_mutsites_overlap_with_iclip_joined_plots(
        Mut_path1, Mut_path2, filt_by,
        Coverage_thresholds, Diff_thresholds,
        Colors, y_lim = c(), y_tick_interval = 500,
        Peak_caller, fig_name = fig_file,
        fig_sz = c(6, 7), sub_plots = sub_plots
      )
    }
  }
}

###############################################################################################################################
# 📝 Suggested Figure Caption
# 
# Figure X. Bar plots showing the number of genes with mutation sites overlapping iCLIP peaks called by five peak callers across UV treatment conditions. Overlap counts were separated by treatment status and plotted for coverage ≥50 and mutation threshold ≥0.03. Dark bars represent genes with overlapping iCLIP peaks; light bars represent the remainder.

