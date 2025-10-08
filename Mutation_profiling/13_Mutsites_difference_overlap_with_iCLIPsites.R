################################################################################
# scatter_plot_utils.R
# Author: Mohamed B
# Description:
# Generates summary scatter plots comparing overlap of mutation sites and regions
# with iCLIP peaks across coverage/difference thresholds and sample pairs.
# Produces publication-ready multi-panel figures.
################################################################################

### 🗂️ Workspace Initialization ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📦 Load Required Libraries ###
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

################################################################################
# 📥 Save Mut-Site/iCLIP Overlap Summary Matrix to .tsv File
################################################################################
write_overlapping_sites <- function(List, List_names, Coverage_thresholds, Diff_thresholds, tsv_file) {
  for (z in seq_along(List)) {
    Mat <- List[[z]]
    Name <- List_names[z]
    
    colnames(Mat) <- Coverage_thresholds
    rownames(Mat) <- Diff_thresholds
    
    df_sub <- melt(Mat, varnames = c("diff_thesh", "cov_thresh"), value.name = Name)
    df <- if (z == 1) df_sub else merge(df, df_sub, by = c("diff_thesh", "cov_thresh"))
  }
  
  write.table(df, file = tsv_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

################################################################################
# 📊 Create Multi-Panel Scatter Plot Comparing Overlaps Between Two Samples
################################################################################
merge_scatter_plots <- function(fig_size, fig_file, Sample1, Sample2, DF,
                                Cov_thresh_colors, diff_thresh_shapes, title) {
  
  ## Columns to use for plotting
  col_select <- list(
    c(paste0(Sample1, "_Sites_overlap_percentage"), paste0(Sample2, "_Sites_overlap_percentage")),
    "sites_diff_overlap", c("sites_diff_overlap", paste0(Sample1, "_Sites")),
    c(paste0(Sample1, "_Regions_overlap_percentage"), paste0(Sample2, "_Regions_overlap_percentage")),
    "regions_diff_overlap", c("regions_diff_overlap", paste0(Sample1, "_Regions"))
  )
  
  ## Axis labels
  x_labs <- c(
    "Difference threshold", "Difference threshold", "Diff in mut sites overlap (%) with iCLIP sites",
    "Difference threshold", "Difference threshold", "Diff in mut regions overlap (%) with iCLIP sites"
  )
  
  y_labs <- c(
    "Mut sites overlap (%) with iCLIP sites",
    "Diff in mut sites overlap (%) with iCLIP sites",
    "Number of Sites",
    "Mut region overlap (%) with iCLIP sites",
    "Diff in mut regions overlap (%) with iCLIP sites",
    "Number of regions"
  )
  
  ## Legend labels
  colour_legends <- rep("Coverage threshold", 6)
  shape_legends <- c("Sample", "", "Difference threshold", "Sample", "", "Difference threshold")
  
  ## Plot config
  Point_sz <- 4
  Alpha <- 0.8
  fig_matrix <- c(2, 3)
  sub_fig_size <- 1 / fig_matrix
  
  fig <- ggdraw()
  ii <- 1
  
  for (i in seq_len(fig_matrix[1])) {
    for (j in seq_len(fig_matrix[2])) {
      if (ii > length(col_select)) break
      
      ## Retrieve column config
      sel <- col_select[[ii]]
      x_lab <- x_labs[ii]
      y_lab <- y_labs[ii]
      color_legend <- colour_legends[ii]
      shape_legend <- shape_legends[ii]
      
      ## Plot selection logic
      if (j %% 3 == 1) {
        # Compare sample1 vs sample2 (shape = sample)
        df1 <- DF[, c("diff_thesh", "cov_thresh", sel[1])]
        df2 <- DF[, c("diff_thesh", "cov_thresh", sel[2])]
        colnames(df1)[3] <- "Z"
        colnames(df2)[3] <- "Z"
        df1$Sample <- Sample1
        df2$Sample <- Sample2
        df <- rbind(df1, df2)
        colnames(df) <- c("x", "Color", "y", "Shape")
        df$Shape <- factor(df$Shape, levels = c(Sample1, Sample2))
        
        p <- plot_scatter3(df, Cov_thresh_colors, Point_sz, Alpha, title, x_lab, y_lab, color_legend, shape_legend)
        
      } else if (j %% 3 == 2) {
        # Standard scatter plot
        df <- DF[, c("diff_thesh", "cov_thresh", sel)]
        colnames(df)[3] <- "Z"
        p <- plot_scatter(df, Cov_thresh_colors, Point_sz, Alpha, title, x_lab, y_lab, color_legend)
        
      } else if (j %% 3 == 0) {
        # Plot with diff_thresh shape mappings
        df <- DF[, c("diff_thesh", "cov_thresh", sel)]
        colnames(df) <- c("Shape", "Color", "x", "y")
        p <- plot_scatter2(df, Cov_thresh_colors, Point_sz, Alpha,
                           diff_thresh_shapes, title, x_lab, y_lab,
                           color_legend, shape_legend)
      }
      
      ## Position plot in layout grid (top-to-bottom)
      fig <- fig + draw_plot(p,
                             x = (i - 1) * sub_fig_size[1],
                             y = (fig_matrix[2] - j) * sub_fig_size[2],
                             width = sub_fig_size[1],
                             height = sub_fig_size[2])
      ii <- ii + 1
    }
  }
  
  ## Save composite figure
  ggsave(filename = fig_file,
         plot = fig,
         width = fig_matrix[1] * fig_size[1],
         height = fig_matrix[2] * fig_size[2])
  
  return(fig)
}

### 📈 1. Plot Linear Scatter Plot — Mutation Overlap by Coverage Threshold ###
plot_scatter <- function(df, Cov_thresh_colors, Point_sz, Alpha,
                         title, x_lab, y_lab, legend_title) {
  
  # Axis padding calculations
  x_range <- range(df$diff_thesh, na.rm = TRUE)
  y_range <- range(df$Z, na.rm = TRUE)
  
  x_min <- x_range[1] - 0.05 * diff(x_range)
  x_max <- x_range[2] + 0.05 * diff(x_range)
  y_min <- floor(min(df$Z) / 5) * 5
  y_max <- ceiling(max(df$Z) / 5) * 5
  x_ticks <- seq(x_min, x_max, by = 0.01)
  y_ticks <- seq(y_min, y_max, by = 5)
  
  # Aesthetic config
  text_size <- 16
  fig <- ggplot(df, aes(x = diff_thesh, y = Z, color = as.factor(cov_thresh))) +
    geom_point(size = Point_sz, alpha = Alpha, stroke = 2) +
    geom_line(linewidth = 1, alpha = Alpha) +
    scale_color_manual(values = Cov_thresh_colors) +
    labs(title = title, x = x_lab, y = y_lab, color = legend_title) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, size = text_size),
      axis.text  = element_text(size = text_size, color = "black"),
      axis.title = element_text(size = text_size),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = "right",
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size)
    ) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_ticks,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = y_ticks,
                       expand = expansion(mult = c(0, 0.02)))
  
  return(fig)
}

### 📈 2. Plot Log-Scaled Scatter — Mutation Rates with Threshold Shape Mapping ###
plot_scatter2 <- function(df, Cov_thresh_colors, Point_sz, Alpha,
                          diff_thresh_shapes, title, x_lab, y_lab,
                          colour_legend, shape_legend) {
  
  # Axis padding & range setup
  x_range <- range(df$x, na.rm = TRUE)
  y_range <- range(df$y, na.rm = TRUE)
  
  x_min <- floor(x_range[1] / 5) * 5
  x_max <- ceiling(x_range[2] / 5) * 5
  x_ticks <- seq(x_min, x_max, by = 5)
  
  y_min <- 10^(floor(log10(min(df$y))))
  y_max <- 10^(ceiling(log10(max(df$y))))
  
  # Aesthetic config
  text_size <- 16
  fig <- ggplot(df, aes(x = x, y = y,
                        color = as.factor(Color),
                        shape = as.factor(Shape))) +
    geom_point(size = Point_sz, alpha = Alpha, stroke = 2) +
    scale_color_manual(values = Cov_thresh_colors) +
    scale_shape_manual(values = diff_thresh_shapes) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "black", linewidth = 0.8) +
    labs(title = title, x = x_lab, y = y_lab,
         color = colour_legend, shape = shape_legend) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, size = text_size),
      axis.text  = element_text(size = text_size, color = "black"),
      axis.title = element_text(size = text_size),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = "right",
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size)
    ) +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_ticks,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_log10(limits = c(y_min, y_max),
                  breaks = 10^seq(0, 7, by = 1),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  expand = expansion(mult = c(0, 0.02)))
  
  return(fig)
}


### Annotated Scatter Plot with Continuous Axes

plot_scatter3 <- function(df, Cov_thresh_colors, Point_sz, Alpha,
                          title, x_lab, y_lab, colour_legend, shape_legend) {
  
  ### ▪ Axis ranges & padding
  x_range <- range(df$x, na.rm = TRUE)
  y_range <- range(df$y, na.rm = TRUE)
  
  x_min <- x_range[1] - 0.05 * diff(x_range)
  x_max <- x_range[2] + 0.05 * diff(x_range)
  y_min <- y_range[1] - 0.05 * diff(y_range)
  y_max <- y_range[2] + 0.05 * diff(y_range)
  
  ### ▪ Aesthetic styling
  text_size <- 16
  Stroke <- 2
  Shapes <- c(16, 17)  # Solid circle, solid triangle
  
  fig <- ggplot(df, aes(x = x, y = y,
                        color = as.factor(Color),
                        shape = as.factor(Shape))) +
    geom_point(size = Point_sz, alpha = Alpha, stroke = Stroke) +
    geom_line(linewidth = 1, alpha = Alpha) +
    scale_color_manual(values = Cov_thresh_colors) +
    scale_shape_manual(values = Shapes) +
    labs(
      title = title,
      x = x_lab, y = y_lab,
      color = colour_legend,
      shape = shape_legend
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.8),
      plot.title       = element_text(hjust = 0.5, size = text_size),
      axis.text        = element_text(size = text_size, color = "black"),
      axis.title       = element_text(size = text_size),
      axis.ticks       = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position  = "right",
      legend.title     = element_text(size = text_size),
      legend.text      = element_text(size = text_size)
    )
  
  ### ▪ Tick interval setup
  x_ticks <- seq(min(df$x), max(df$x), by = 0.01)
  y_ticks <- seq(floor(min(df$y) / 5) * 5, ceiling(max(df$y) / 5) * 5, by = 5)
  
  fig <- fig +
    scale_x_continuous(limits = c(x_min, x_max), breaks = x_ticks,
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = y_ticks,
                       expand = expansion(mult = c(0, 0.02)))
  
  return(fig)
}


################################################################################

### ▪ Setup: Thresholds, Styles, Sample Paths ###
Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)

Cov_thresh_colors <- c("#000000", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
names(Cov_thresh_colors) <- Coverage_thresholds

diff_thresh_shapes <- c(16, 1, 15, 0, 17, 2)
names(diff_thresh_shapes) <- Diff_thresholds

Samples1 <- "SRSF3"
Samples2 <- "P19_WT_no_UV"

filt_by <- "FDR"  # Can switch to "pvalue" if needed

Peak_callers <- c("CLIPper", "iCount", "pureclip", "piranha", "paraclu")
iclip_paths <- c(
  "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/clipper/SRSF3_merged.peaks.bed.gz",
  "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/icount/SRSF3_merged.3nt_3nt.peaks.bed.gz",
  "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/pureclip/SRSF3_merged.8nt.peaks.bed.gz",
  "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/piranha/SRSF3_merged.3nt_3nt.peaks.bed.gz",
  "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/paraclu/SRSF3_merged.10_200nt_2.peaks.bed.gz"
)

paths1 <- rep("Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/", 2)
ref_subfolder1 <- c("ref_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV")

paths2 <- rep("Results/Nanopore/4_noUV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/P19_WT_TotalRNA_cDNA_no_UV_merge/", 2)
ref_subfolder2 <- c("ref_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV")

subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

### ▪ Main Loop ###
for (ii in seq_along(Peak_callers)) {
  Peak_caller <- Peak_callers[ii]
  iclip_path <- iclip_paths[ii]
  
  icLIP_GR <- load_CLIP_tags(iclip_path, Peak_caller, Left_offset = 0, Right_offset = 0)
  
  for (subfolder1 in subfolders1) {
    for (subfolder2 in subfolders2) {
      Mut_paths1 <- paste0(paths1, subfolder1, "/", subfolder2, "/", ref_subfolder1, "/", filt_by, "/")
      Mut_paths2 <- paste0(paths2, subfolder1, "/", subfolder2, "/", ref_subfolder2, "/", filt_by, "/")
      
      for (i in seq_along(Mut_paths1)) {
        Mut_path1 <- Mut_paths1[i]
        Mut_path2 <- Mut_paths2[i]
        tsv_file <- paste0(Peak_caller, "_iCLIPsites_overlap_with_Pos_Mutsites.tsv")
        
        ### ▪ Load overlap data
        Olap1 <- read.delim(file.path(Mut_path1, tsv_file))
        Olap2 <- read.delim(file.path(Mut_path2, tsv_file))
        
        colnames(Olap1)[-(1:2)] <- paste0(Samples1, "_", colnames(Olap1)[-(1:2)])
        colnames(Olap2)[-(1:2)] <- paste0(Samples2, "_", colnames(Olap2)[-(1:2)])
        
        DF <- merge(Olap1, Olap2, by = c("diff_thesh", "cov_thresh"))
        
        ### ▪ Compute differences
        DF$sites_diff_overlap <- DF[[paste0(Samples1, "_Sites_overlap_percentage")]] - DF[[paste0(Samples2, "_Sites_overlap_percentage")]]
        DF$regions_diff_overlap <- DF[[paste0(Samples1, "_Regions_overlap_percentage")]] - DF[[paste0(Samples2, "_Regions_overlap_percentage")]]
        
        ### ▪ Initialize matrices
        matrix_template <- matrix(NA, nrow = length(Diff_thresholds), ncol = length(Coverage_thresholds),
                                  dimnames = list(paste0("T_", Diff_thresholds), paste0("T_", Coverage_thresholds)))
        
        Sites_overlap <- matrix_template
        sites_diff_overlap <- matrix_template
        regions_overlap <- matrix_template
        regions_diff_overlap <- matrix_template
        
        ### ▪ Populate matrices
        for (k in seq_along(Coverage_thresholds)) {
          for (m in seq_along(Diff_thresholds)) {
            id <- DF$diff_thesh == Diff_thresholds[m] & DF$cov_thresh == Coverage_thresholds[k]
            if (any(id)) {
              Sites_overlap[m, k] <- DF[[paste0(Samples1, "_Sites_overlap_percentage")]][id]
              sites_diff_overlap[m, k] <- DF$sites_diff_overlap[id]
              regions_overlap[m, k] <- DF[[paste0(Samples1, "_Regions_overlap_percentage")]][id]
              regions_diff_overlap[m, k] <- DF$regions_diff_overlap[id]
            }
          }
        }
        
        ### ▪ Generate scatter plot
        fig_file <- file.path(Mut_path1, paste0("Mutsites_diff_in_overlap_with_", Peak_caller, "_iCLIPsites_scatter.pdf"))
        f2 <- merge_scatter_plots(
          fig_size = c(7.5, 7),
          fig_file = fig_file,
          Sample1 = Samples1, Sample2 = Samples2,
          DF = DF,
          Cov_thresh_colors = Cov_thresh_colors,
          diff_thresh_shapes = diff_thresh_shapes,
          title = Peak_caller
        )
      }
    }
  }
}

#######################################################################################
# 📝 Supplementary Figure Caption Suggestion
# 
# Multi-panel scatter plots showing mutation-iCLIP site/region overlap percentages across varying coverage and difference thresholds. Comparative panels illustrate sample-specific changes, differential overlap scores, and absolute mutation counts.

# 📝 Figure Caption Suggestions
# 
# Fig. X. Scatter plots showing overlap of mutation sites and regions with iCLIP peaks across coverage and difference thresholds. Panels display either absolute counts or mutation rates (log-scaled), with visual distinction by threshold and sample condition.

# 🖼️ Figure Caption Suggestion
# 
# Fig. X. Comparative scatter plots of mutation site overlaps with iCLIP peaks called by various algorithms. Each panel shows pairwise differences between SRSF3 and P19_WT_no_UV across coverage and difference thresholds. Colors encode coverage levels, shapes represent difference thresholds.


