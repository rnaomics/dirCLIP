################################################################################
# scatter_plot_analysis.R
# Author: Mohamed B
# Description:
# Generates comparative scatter plots for mutation site/region overlaps
# across multiple coverage and difference thresholds grouped by conditions.
################################################################################

### 📁 Setup
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📚 Libraries
library(Rsubread)
library(writexl)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggcorrplot)
library(scales)
library(GenomicRanges)
library(S4Vectors)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(cowplot)

######################################################################################
### 🔹 Plot Individual Scatter Plot
plot_scatter1 <- function(df, Cond_colors, Point_sz, Alpha,
                          diff_thresh_shapes, title, x_lab, y_lab,
                          x_tick_interval, colour_legend, shape_legend) {
  
  x_min <- min(df$x, na.rm = TRUE)
  x_max <- max(df$x, na.rm = TRUE)
  x_range <- x_max - x_min
  x_min <- x_min - x_range * 0.05
  x_max <- x_max + x_range * 0.05
  
  y_min <- min(df$y, na.rm = TRUE)
  y_max <- max(df$y, na.rm = TRUE)
  y_range <- y_max - y_min
  y_min <- y_min * 0.8
  y_max <- y_max + y_range * 0.05
  
  text_size <- 16
  
  fig <- ggplot(df, aes(x = x, y = y, color = as.factor(Color))) +
    geom_point(size = Point_sz, alpha = Alpha) +
    geom_line(linewidth = 1, alpha = Alpha) +
    scale_color_manual(values = Cond_colors) +
    labs(title = title, x = x_lab, y = y_lab, color = colour_legend) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      axis.line  = element_line(color = "black", linewidth = 0.8),
      axis.text  = element_text(size = text_size, color = "black"),
      axis.title = element_text(size = text_size),
      plot.title = element_text(hjust = 0.5, size = text_size),
      legend.position = "right",
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = seq(floor(x_min / x_tick_interval) * x_tick_interval,
                   ceiling(x_max / x_tick_interval) * x_tick_interval,
                   by = x_tick_interval),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_log10(
      breaks = 10^seq(0, 7, by = 1),
      labels = trans_format("log10", math_format(10^.x)),
      expand = expansion(mult = c(0, 0.02))
    )
  
  return(fig)
}

### 🔸 Merge All Panels into Composite Plot
merge_scatter_plots <- function(fig_size, fig_file, Sample, DF, Condition_colors) {
  
  cols   <- c("Sites", "Sites", "Regions", "Regions")
  thresh <- c("cov_thresh", "diff_thresh", "cov_thresh", "diff_thresh")
  vals   <- c(50, 0.03, 50, 0.03)
  
  x_ticks <- c(0.01, 100, 0.01, 100)
  x_labs  <- c("Difference threshold", "Coverage threshold",
               "Difference threshold", "Coverage threshold")
  y_labs  <- c("Number of sites", "Number of sites",
               "Number of regions", "Number of regions")
  
  fig_matrix <- c(2, 2)
  sub_sz <- 1 / fig_matrix
  fig <- ggdraw()
  ii <- 1
  
  for (i in 1:fig_matrix[1]) {
    for (j in 1:fig_matrix[2]) {
      df_sub <- if (thresh[ii] == "cov_thresh") {
        df <- DF[DF$cov_thresh == vals[ii], ]
        df[, c("diff_thesh", "Condition", cols[ii])]
      } else {
        df <- DF[DF$diff_thesh == vals[ii], ]
        df[, c("cov_thresh", "Condition", cols[ii])]
      }
      
      colnames(df_sub) <- c("x", "Color", "y")
      title <- paste0(Sample, " : ", thresh[ii], " = ", vals[ii])
      
      p <- plot_scatter1(df_sub, Condition_colors, Point_sz = 4, Alpha = 1,
                         diff_thresh_shapes = NULL,
                         title, x_labs[ii], y_labs[ii], x_ticks[ii],
                         colour_legend = "Conditions", shape_legend = "")
      
      fig <- fig + draw_plot(p,
                             x = (i - 1) * sub_sz[1],
                             y = (fig_matrix[2] - j) * sub_sz[2],
                             width  = sub_sz[1],
                             height = sub_sz[2])
      ii <- ii + 1
    }
  }
  
  ggsave(fig_file, plot = fig,
         width  = fig_size[1] * fig_matrix[1],
         height = fig_size[2] * fig_matrix[2])
  return(fig)
}

###############################################################################################################################

### 🔁 Input Paths and Metadata
Paths <- c(
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_TotalRNA_cDNA_no_UV/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/1_No_del/ref_TotalRNA_cDNA_no_UV/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_TotalRNA_cDNA_UV/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/1_No_del/ref_TotalRNA_cDNA_UV/"
)

Conditions <- c("-UV, del", "-UV, nodel", "+UV, del", "+UV, nodel")
names(Conditions) <- Conditions
Sample <- "SRSF3"
Result_path = "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_TotalRNA_cDNA_no_UV/"

filt_by <- "FDR"

### 🔎 iCLIP Metadata
Peak_caller <- "iCount"
tsv_file <- paste0(Peak_caller, "_iCLIPsites_overlap_with_Pos_Mutsites.tsv")

### 🎨 Condition Colors
Condition_colors <- c("#333333", "#666666", "#999999", "#CCCCCC")
names(Condition_colors) <- Conditions

### 🧮 Load and Combine Data
DF <- do.call(rbind, lapply(seq_along(Paths), function(i) {
  Mut_path <- file.path(Paths[i], filt_by)
  df <- read.delim(file.path(Mut_path, tsv_file))
  df$Condition <- Conditions[i]
  return(df)
}))

### ✨ Generate Figure
merge_scatter_plots(fig_size = c(7.5, 7),
                    fig_file = file.path(Result_path, filt_by, "Mutsites_vs_cov_diff_thresh.pdf"),
                    Sample = Sample, DF = DF, Condition_colors = Condition_colors)

##############################################################################################################

