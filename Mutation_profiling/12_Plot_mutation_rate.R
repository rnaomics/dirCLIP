################################################################################
# mutation_rate_boxplot.R
# Author: Mohamed B
# Description:
# Processes mutation rate files from Nanopore DIRCLIP, filters by coverage
# and iCLIP overlap, calculates boxplot statistics, and generates custom-styled
# box plots for comparative visualization across samples and deletion status.
################################################################################

### 💼 Workspace Setup ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📦 Load Required Libraries ###
library(Rsubread)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(S4Vectors)
library(matrixStats)
library(GenomicRanges)

################################################################################
# 🔍 1. Filter Mutation Rate Table by iCLIP Peak Overlap
################################################################################
get_overlap_with_iCLIP <- function(Mut_Rate, iCLIP_GR) {
  Mut_Rate_GR <- GRanges(seqnames = Mut_Rate$Chr,
                         ranges = IRanges(start = Mut_Rate$loc, end = Mut_Rate$loc),
                         strand = "*")
  overlaps <- findOverlaps(Mut_Rate_GR, iCLIP_GR)
  Mut_Rate_filtered <- Mut_Rate[unique(queryHits(overlaps)), ]
  return(Mut_Rate_filtered)
}

################################################################################
# 📊 2. Calculate Boxplot Statistics Manually
################################################################################
get_boxplot_values <- function(Values) {
  Q1 <- quantile(Values, 0.25)
  Q2 <- quantile(Values, 0.5)
  Q3 <- quantile(Values, 0.75)
  IQR <- Q3 - Q1
  
  lower_whisker <- min(Values[Values >= (Q1 - 1.5 * IQR)])
  upper_whisker <- max(Values[Values <= (Q3 + 1.5 * IQR)])
  N <- length(Values)
  
  return(data.frame(lower_whisker, Q1, Q2, Q3, upper_whisker, N))
}

################################################################################
# 📈 3. Boxplot Stats Wrapper Over Coverage Thresholds
################################################################################
get_boxplot_values_for_Coverage_thresholds <- function(Paths, Samples, isdel,
                                                       Mut_file, Coverage_thresholds,
                                                       Common_locs, iCLIP_GR) {
  box_plot_values_list <- list()
  
  for (j in seq_along(Paths)) {
    Mut_Rate <- read.delim(file.path(Paths[j], Mut_file), header = FALSE)
    colnames(Mut_Rate) <- c("Chr", "loc", "Coverage", "Mut_rate")
    
    if (length(Common_locs) > 0) {
      Mut_Rate <- merge(Mut_Rate, Common_locs, by = c("Chr", "loc"))
    }
    
    if (is(iCLIP_GR, "GRanges")) {
      Mut_Rate <- get_overlap_with_iCLIP(Mut_Rate, iCLIP_GR)
    }
    
    box_plot_values <- do.call(rbind, lapply(Coverage_thresholds, function(cov_thresh) {
      Values <- Mut_Rate[Mut_Rate$Coverage >= cov_thresh, "Mut_rate"]
      get_boxplot_values(Values)
    }))
    
    rownames(box_plot_values) <- Coverage_thresholds
    box_plot_values_list[[j]] <- box_plot_values
  }
  
  return(tibble(box_plot_values_list, Samples, isdel))
}

################################################################################
# 📊 4. Plot Styled Boxplots (Grouped by Sample)
################################################################################
make_box_plot_type1 <- function(bp_values, Colors, x_lab, y_lab, Title,
                                y_limit, text_size, text_size2) {
  fig <- ggplot(bp_values, aes(x = Samples, fill = Samples)) +
    geom_boxplot(aes(ymin = lower_whisker, lower = Q1, middle = Q2, upper = Q3, ymax = upper_whisker),
                 stat = "identity", width = 0.6, color = "black", linewidth = 0.6, alpha = 0.8) +
    scale_fill_manual(values = Colors) +
    geom_text(aes(x = Samples, y = upper_whisker * 1.05, label = paste0("n=", N)),
              vjust = 0, size = text_size2, fontface = "bold", color = "black") +
    labs(y = y_lab, x = x_lab, fill = x_lab, title = Title) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.8),
      axis.text.x       = element_text(size = text_size, color = "black", angle = 45, hjust = 1),
      axis.text.y       = element_text(size = text_size, color = "black"),
      axis.title        = element_text(size = text_size),
      legend.position   = "none",
      legend.title      = element_text(size = text_size),
      legend.text       = element_text(size = text_size),
      axis.ticks        = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    )
  return(fig)
}

################################################################################
# 📊 5. Plot Styled Boxplots (Grouped by Deletion Status)
################################################################################
make_box_plot_type2 <- function(bp_values, Colors, x_lab, y_lab, Title,
                                legend_title, y_limit, text_size, text_size2) {
  fig <- ggplot(bp_values, aes(x = is_del, fill = Samples)) +
    geom_boxplot(aes(ymin = lower_whisker, lower = Q1, middle = Q2, upper = Q3, ymax = upper_whisker),
                 stat = "identity", width = 0.6, color = "black", linewidth = 0.6, alpha = 0.8) +
    scale_fill_manual(values = Colors) +
    labs(y = y_lab, x = x_lab, fill = legend_title, title = Title) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.8),
      axis.text.x       = element_text(size = text_size, color = "black"),
      axis.text.y       = element_text(size = text_size, color = "black"),
      axis.title        = element_text(size = text_size),
      legend.title      = element_text(size = text_size),
      legend.text       = element_text(size = text_size),
      axis.ticks        = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    )
  return(fig)
}

################################################################################
# Generate grouped box plots of mutation rates from Nanopore datasets.
# Filters samples by deletion status and overlays coverage threshold comparisons.
################################################################################

### 📊 Function to Create Box Plot for Specific Coverage Threshold ###
make_box_plot <- function(box_plot_values_list, cov_thresh, Samples, isdel,
                          Sample_order, Sample_Colors, del_Colors,
                          x_lab, y_lab, Title, y_lim, y_tick_interval) {
  
  selected <- box_plot_values_list[box_plot_values_list$isdel %in% isdel, ]
  Samples <- selected$Samples
  is_del  <- selected$isdel
  
  bp_values <- do.call(rbind, lapply(seq_along(Samples), function(i) {
    selected[[i]]$box_plot_values_list[[1]][as.character(cov_thresh), ]
  }))
  
  bp_values$Samples <- factor(Samples, levels = Sample_order)
  bp_values$is_del  <- is_del
  
  ## Set default axis limits
  if (length(y_lim) == 0) {
    y_lim <- c(0, max(bp_values$upper_whisker) * 1.1)
  }
  
  ## Determine which plot style to use
  text_size <- 16
  text_size2 <- 6
  
  if (length(unique(is_del)) == 1) {
    fig <- make_box_plot_type1(bp_values, Sample_Colors, x_lab, y_lab, Title, y_lim, text_size, text_size2)
  } else {
    fig <- make_box_plot_type2(bp_values, Sample_Colors, x_lab, y_lab, Title, legend_title = "Sample", y_lim, text_size, text_size2)
  }
  
  ## Configure Y-axis ticks
  y_min <- floor(y_lim[1] / y_tick_interval) * y_tick_interval
  y_max <- ceiling(y_lim[2] / y_tick_interval) * y_tick_interval
  y_ticks <- seq(y_min, y_max, by = y_tick_interval)
  
  fig <- fig + scale_y_continuous(limits = c(y_min, y_max), breaks = y_ticks,
                                  expand = expansion(mult = c(0, 0.02)))
  
  return(fig)
}

### 🎯 Function to Merge Multiple Box Plots Across Coverage Thresholds ###
make_joined_box_plots <- function(box_plot_values_list, Coverage_thresholds,
                                  Samples, isdel, Sample_order,
                                  Sample_Colors, del_Colors,
                                  x_lab, y_lab, y_lim, y_tick_interval,
                                  fig_name, fig_sz, sub_plots) {
  
  fig <- ggdraw()
  sub_fig_sz <- 1 / sub_plots
  grid_pos <- expand.grid(x = 0:(sub_plots[1]-1), y = 0:(sub_plots[2]-1))
  
  for (i in seq_along(Coverage_thresholds)) {
    cov_thresh <- Coverage_thresholds[i]
    print(cov_thresh)
    
    f <- make_box_plot(box_plot_values_list, cov_thresh, Samples, isdel,
                       Sample_order, Sample_Colors, del_Colors,
                       x_lab, y_lab, paste0("Cov thresh = ", cov_thresh),
                       y_lim, y_tick_interval)
    
    fig <- fig + draw_plot(f,
                           x = grid_pos$x[i] * sub_fig_sz[1],
                           y = grid_pos$y[i] * sub_fig_sz[2],
                           width = sub_fig_sz[1],
                           height = sub_fig_sz[2])
  }
  
  ## Save PDF
  ggsave(filename = fig_name, plot = fig,
         width = fig_sz[1] * sub_plots[1],
         height = fig_sz[2] * sub_plots[2],
         limitsize = FALSE, device = cairo_pdf)
}

#############################################################################
# Define file paths and sample conditions
Paths <- c(
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/subsample1/1_No_del/",
  "Results/Nanopore/4_noUV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/TotalRNA_cDNA_no_UV_merge/q_8_Q_8/subsample1/1_No_del/",
  "Results/Nanopore/5_UV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/TotalRNA_cDNA_UV_merge/q_8_Q_8/subsample1/1_No_del/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/subsample1/2_Inc_del/",
  "Results/Nanopore/4_noUV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/TotalRNA_cDNA_no_UV_merge/q_8_Q_8/subsample1/2_Inc_del/",
  "Results/Nanopore/5_UV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/TotalRNA_cDNA_UV_merge/q_8_Q_8/subsample1/2_Inc_del/"
)

Samples <- c("SRSF3 CLIP", "Total RNA non-UV", "Total RNA UV",
             "SRSF3 CLIP", "Total RNA non-UV", "Total RNA UV")

Sample_order <- c("Total RNA non-UV", "Total RNA UV", "SRSF3 CLIP")
isdel <- c("No", "No", "No", "Yes", "Yes", "Yes")

Mut_file <- "Mutation_rate_cov_10_mut_0.tsv.gz"

Sample_Colors <- c("SRSF3 CLIP" = "#66c2a5",
                   "Total RNA non-UV" = "#08519C",
                   "Total RNA UV" = "#9ECAE1")

del_Colors <- c("Yes" = "gray40", "No" = "gray70")

# Output paths
Result_path_nodel = "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/1_No_del/"
Result_path_del = "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/"

Coverage_thresholds <- c(50)  # You can add more thresholds here

# Get box plot data
box_plot_values_list1 <- get_boxplot_values_for_Coverage_thresholds(
  Paths, Samples, isdel, Mut_file, Coverage_thresholds, Common_locs = c(), iCLIP_GR = "")

# Save plots (no deletions)
fig_sz <- c(5, 7)
y_tick_interval <- 0.01
dir.create(Result_path_nodel, recursive = TRUE, showWarnings = FALSE)

make_joined_box_plots(box_plot_values_list1, Coverage_thresholds, Samples,
                      isdel = "No", Sample_order, Sample_Colors, del_Colors,
                      x_lab = "Sample", y_lab = "Mut rate",
                      y_lim = c(), y_tick_interval,
                      fig_name = file.path(Result_path_nodel, "Mut_Rate_nodel.pdf"),
                      fig_sz = fig_sz, sub_plots = c(1,1))

# Save plots (with deletions)
y_tick_interval <- 0.025
dir.create(Result_path_del, recursive = TRUE, showWarnings = FALSE)

make_joined_box_plots(box_plot_values_list1, Coverage_thresholds, Samples,
                      isdel = "Yes", Sample_order, Sample_Colors, del_Colors,
                      x_lab = "Sample", y_lab = "Mut rate",
                      y_lim = c(), y_tick_interval,
                      fig_name = file.path(Result_path_del, "Mut_Rate_del.pdf"),
                      fig_sz = fig_sz, sub_plots = c(1,1))

# Save combined plots
fig_sz <- c(7, 7)
make_joined_box_plots(box_plot_values_list1, Coverage_thresholds, Samples,
                      isdel = c("No", "Yes"), Sample_order, Sample_Colors, del_Colors,
                      x_lab = "Dels in mutation", y_lab = "Mut rate",
                      y_lim = c(), y_tick_interval,
                      fig_name = file.path(Result_path_del, "Mut_Rate_del_nodel.pdf"),
                      fig_sz = fig_sz, sub_plots = c(1,1))

################################################################################################################
# 💡 Suggested Caption for Supplementary Figure
# 
# Boxplot summary of mutation rates for positions overlapping with iCLIP peaks. Mutation sites were filtered by coverage thresholds and grouped either by sample identity or deletion status. Whiskers calculated using 1.5×IQR rule.




