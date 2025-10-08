################################################################################
# mutation_overlap_analysis.R
# Author: Mohamed B
# Description: End-to-end pipeline for assessing mutation site overlaps with
# iCLIP data and visualizing results.
################################################################################

# Set working directory and library path
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

# Load necessary libraries
library(Rsubread)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(GenomicRanges)
library(S4Vectors)
library(reshape2)
library(RColorBrewer)

################################################################################
# Function Definitions
################################################################################

# Load iCLIP tag data and convert to GRanges object
load_CLIP_tags <- function(iclip_path, Peak_caller, Left_offset, Right_offset){
  if (Peak_caller == "Michael_unknown"){
    iclip_bed <- read.table(iclip_path, skip = 1)
    iclip_bed$Strand <- ifelse(iclip_bed$V4 < 0, "-", "+")
    iclip_bed$V2[iclip_bed$Strand == "+"] <- iclip_bed$V2[iclip_bed$Strand == "+"] - Left_offset
    iclip_bed$V3[iclip_bed$Strand == "+"] <- iclip_bed$V3[iclip_bed$Strand == "+"] + Right_offset
    iclip_bed$V2[iclip_bed$Strand == "-"] <- iclip_bed$V2[iclip_bed$Strand == "-"] - Right_offset
    iclip_bed$V3[iclip_bed$Strand == "-"] <- iclip_bed$V3[iclip_bed$Strand == "-"] + Left_offset
    iclip_GR <- GRanges(seqnames = iclip_bed$V1, ranges = IRanges(iclip_bed$V2, iclip_bed$V3), strand = iclip_bed$Strand)
  } else {
    iclip_bed <- read.delim(iclip_path, header = FALSE)
    iclip_bed$V2[iclip_bed$V6 == "+"] <- iclip_bed$V2[iclip_bed$V6 == "+"] - Left_offset
    iclip_bed$V3[iclip_bed$V6 == "+"] <- iclip_bed$V3[iclip_bed$V6 == "+"] + Right_offset
    iclip_bed$V2[iclip_bed$V6 == "-"] <- iclip_bed$V2[iclip_bed$V6 == "-"] - Right_offset
    iclip_bed$V3[iclip_bed$V6 == "-"] <- iclip_bed$V3[iclip_bed$V6 == "-"] + Left_offset
    iclip_GR <- GRanges(seqnames = iclip_bed$V1, ranges = IRanges(iclip_bed$V2, iclip_bed$V3), strand = iclip_bed$V6)
  }
  return(iclip_GR)
}

# Load mutation sites and convert to GRanges
load_sites <- function(Path, Left_offset, Right_offset){
  Sites <- read.delim(Path, header = FALSE)
  Sites <- Sites[!is.na(Sites$V2), ]
  Sites <- Sites[!grepl("random", Sites$V1), ]
  Sites$strand <- "*"
  Sites_GR <- GRanges(seqnames = Sites$V1, ranges = IRanges(Sites$V2 - Left_offset, Sites$V2 + Right_offset), strand = Sites$strand)
  return(Sites_GR)
}

# Load mutation regions and convert to GRanges
load_regions <- function(Path, Left_offset, Right_offset){
  Sites <- read.delim(Path, header = TRUE)
  Sites <- Sites[!is.na(Sites$start), ]
  Sites <- Sites[!grepl("random", Sites$chr), ]
  Sites$strand <- "*"
  Sites_GR <- GRanges(seqnames = Sites$chr, ranges = IRanges(Sites$start - Left_offset, Sites$end + Right_offset), strand = Sites$strand)
  return(Sites_GR)
}

# Compute percent overlap between GRanges objects
get_overlap_percentage <- function(GR_A, GR_B){
  id <- findOverlaps(GR_A, GR_B)
  percent <- round(length(unique(queryHits(id))) / length(GR_A) * 100, 2)
  return(percent)
}

# Write overlap results to a TSV file
write_overlapping_sites <- function(List, List_names, Coverage_thresholds, Diff_thresholds, tsv_file){
  for (z in seq_along(List)){
    Mat <- List[[z]]
    Name <- List_names[z]
    colnames(Mat) <- Coverage_thresholds
    rownames(Mat) <- Diff_thresholds
    melted <- melt(Mat, varnames = c("diff_thresh", "cov_thresh"), value.name = Name)
    if (z == 1){
      df <- melted
    } else {
      df <- merge(df, melted, by = c("diff_thresh", "cov_thresh"))
    }
  }
  write.table(df, tsv_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

# Plot matrix of correlation or overlap metrics
plot_matrix <- function(corr, method, range = NULL, lab = TRUE, title = "", x_lab = "", y_lab = "", legend_title = ""){
  corr_df <- melt(corr)
  colnames(corr_df) <- c("Var1", "Var2", "value")
  corr_df$Var1 <- as.factor(corr_df$Var1)
  corr_df$Var2 <- as.factor(corr_df$Var2)
  
  if (is.null(range)) range <- range(corr_df$value)
  
  p <- ggplot(corr_df, aes(x = Var1, y = Var2, fill = value)) +
    theme_minimal(base_size = 14) +
    labs(x = x_lab, y = y_lab, fill = legend_title, title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          panel.grid = element_blank()) +
    scale_fill_gradientn(colors = colorRampPalette(brewer.pal(8, "Blues"))(20), limits = range)
  
  if (method == "circle"){
    p <- p + geom_point(aes(size = value), shape = 21, color = "gray") +
      scale_size(range = c(4, 10)) +
      guides(size = "none")
  } else {
    p <- p + geom_tile(color = "gray")
  }
  
  if (lab){
    label_vals <- round(corr_df$value, digits = 2)
    p <- p + geom_text(aes(label = label_vals), color = "black", size = 4)
  }
  
  return(p)
}

# Combine multiple matrix plots into a single figure
merge_matrix_plots <- function(fig_matrix, fig_size, fig_file, DF, method, range, lab, title, x_lab, y_lab, legend_titles){
  fig <- ggdraw()
  idx <- 1
  for (i in 1:fig_matrix[1]){
    for (j in 1:fig_matrix[2]){
      if (idx <= length(DF)){
        p <- plot_matrix(corr = DF[[idx]], method = method, range = range,
                         lab = lab, title = title[idx], x_lab = x_lab,
                         y_lab = y_lab, legend_title = legend_titles)
        fig <- fig + draw_plot(p, x = (i - 1)/fig_matrix[1], y = (fig_matrix[2] - j)/fig_matrix[2],
                               width = 1/fig_matrix[1], height = 1/fig_matrix[2])
      }
      idx <- idx + 1
    }
  }
  ggsave(fig_file, plot = fig, width = fig_matrix[1]*fig_size[1], height = fig_matrix[2]*fig_size[2])
  return(fig)
}


###########################################################################################################################

# Define input mutation paths and reference folders
paths <- c("Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
           "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
           "Results/Nanopore/4_noUV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/P19_WT_TotalRNA_cDNA_no_UV_merge/")

ref_subfolder <- c("ref_TotalRNA_cDNA_no_UV",
                   "ref_P19_WT_TotalRNA_cDNA_no_UV",
                   "ref_TotalRNA_cDNA_no_UV")

subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

# Define peak callers and corresponding iCLIP BED file paths
Peak_callers <- c("CLIPper", "iCount", "pureclip", "piranha", "paraclu")

iclip_paths <- c("Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/clipper/SRSF3_merged.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/icount/SRSF3_merged.3nt_3nt.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/pureclip/SRSF3_merged.8nt.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/piranha/SRSF3_merged.3nt_3nt.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/paraclu/SRSF3_merged.10_200nt_2.peaks.bed.gz")

# Thresholds for filtering mutation calls
Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)

# Filtering method and suffix for subfolder paths
filt_by <- "FDR"
sub_folder_suffix <- ""

################################################################################
# Main analysis loop across peak callers
################################################################################

for (ii in seq_along(Peak_callers)) {
  Peak_caller <- Peak_callers[ii]
  iclip_path <- iclip_paths[ii]
  
  # Load iCLIP GRanges object
  icLIP_GR <- load_CLIP_tags(iclip_path, Peak_caller, Left_offset = 0, Right_offset = 0)
  
  for (subfolder1 in subfolders1) {
    for (subfolder2 in subfolders2) {
      
      # Assemble mutation profile paths
      Mut_paths <- paste0(paths, subfolder1, "/", subfolder2, "/", ref_subfolder, "/")
      
      for (i in seq_along(Mut_paths)) {
        print(Mut_paths[i])
        
        # Initialize matrices to store overlap metrics
        Pos_sites_Overlap_pecentage <- matrix(NA, nrow = length(Diff_thresholds), ncol = length(Coverage_thresholds),
                                              dimnames = list(Diff_thresholds, Coverage_thresholds))
        
        Pos_regions_Overlap_pecentage <- Pos_sites_Overlap_pecentage
        iCLIP_Pos_sites_overlap_percentage <- Pos_sites_Overlap_pecentage
        iCLIP_Pos_regions_overlap_percentage <- Pos_sites_Overlap_pecentage
        N_Pos_sites <- Pos_sites_Overlap_pecentage
        N_Pos_regions <- Pos_sites_Overlap_pecentage
        
        # Loop through coverage × diff thresholds
        for (k in seq_along(Coverage_thresholds)) {
          for (m in seq_along(Diff_thresholds)) {
            
            cov_thresh <- Coverage_thresholds[k]
            diff_thresh <- Diff_thresholds[m]
            
            Mut_path <- paste0(Mut_paths[i], sub_folder_suffix, "cov_thresh_", cov_thresh, "/", filt_by, "/")
            pos_mut_file <- paste0("Pos_mut_thresh_", diff_thresh, ".tsv")
            pos_mut_region <- paste0("Pos_mut_thresh_", diff_thresh, "_region.tsv")
            
            # Load mutation GRanges
            Pos_Mut_sites_GR <- load_sites(paste0(Mut_path, pos_mut_file), Left_offset = 5, Right_offset = 5)
            Pos_Mut_regions_GR <- load_regions(paste0(Mut_path, pos_mut_region), Left_offset = 5, Right_offset = 5)
            
            # Compute overlaps
            Pos_sites_Overlap_pecentage[m, k] <- get_overlap_percentage(Pos_Mut_sites_GR, icLIP_GR)
            Pos_regions_Overlap_pecentage[m, k] <- get_overlap_percentage(Pos_Mut_regions_GR, icLIP_GR)
            
            iCLIP_Pos_sites_overlap_percentage[m, k] <- get_overlap_percentage(icLIP_GR, Pos_Mut_sites_GR)
            iCLIP_Pos_regions_overlap_percentage[m, k] <- get_overlap_percentage(icLIP_GR, Pos_Mut_regions_GR)
            
            # Store mutation counts
            N_Pos_sites[m, k] <- length(Pos_Mut_sites_GR)
            N_Pos_regions[m, k] <- length(Pos_Mut_regions_GR)
          }
        }
        
        # Collect matrices for plotting
        Pos_Overlap_list <- list(t(N_Pos_sites),
                                 t(Pos_sites_Overlap_pecentage),
                                 t(N_Pos_regions),
                                 t(Pos_regions_Overlap_pecentage))
        
        Pos_overlap_legends <- c("No. of Mut sites",
                                 "Mut sites overlap with iCLIP sites (%)",
                                 "No. of Mut regions",
                                 "Mut regions overlap with iCLIP sites (%)")
        
        iCLIP_Mut_overlap_percentage_list <- list(t(iCLIP_Pos_sites_overlap_percentage),
                                                  t(iCLIP_Pos_regions_overlap_percentage))
        
        iCLIP_Mut_overlap_percentage_legends <- c("iCLIP sites overlap with Mut regions (%)",
                                                  "iCLIP sites overlap with Mut regions (%)")
        
        Mut_Score_list <- list(t(Pos_sites_Overlap_pecentage * iCLIP_Pos_sites_overlap_percentage / 100),
                               t(Pos_regions_Overlap_pecentage * iCLIP_Pos_regions_overlap_percentage / 100))
        
        Mut_Score_legends <- c("Mut sites score",
                               "Mut regions score")
        
        # Create output folder if needed
        Path2 <- paste0(Mut_paths[i], sub_folder_suffix, filt_by, "/")
        if (!dir.exists(Path2)) dir.create(Path2, recursive = TRUE)
        
        # Plot overlap heatmaps
        merge_matrix_plots(c(2, 2), c(7, 5),
                           fig_file = paste0(Path2, "Pos_Mutsites_overlap_with_", Peak_caller, "_iCLIPsites.pdf"),
                           DF = Pos_Overlap_list,
                           method = "square", range = NULL, lab = TRUE,
                           title = paste0(Peak_caller, " : ", Pos_overlap_legends),
                           x_lab = "Coverage threshold", y_lab = "Difference threshold",
                           legend_titles = "")
        
        merge_matrix_plots(c(2, 2), c(7, 5),
                           fig_file = paste0(Path2, Peak_caller, "_iCLIPsites_overlap_with_Pos_Mutsites.pdf"),
                           DF = iCLIP_Mut_overlap_percentage_list,
                           method = "square", range = NULL, lab = TRUE,
                           title = paste0(Peak_caller, " : ", iCLIP_Mut_overlap_percentage_legends),
                           x_lab = "Coverage threshold", y_lab = "Difference threshold",
                           legend_titles = "")
        
        # Save overlap tables
        write_overlapping_sites(List = list(N_Pos_sites,
                                            Pos_sites_Overlap_pecentage,
                                            N_Pos_regions,
                                            Pos_regions_Overlap_pecentage),
                                List_names = c("Sites",
                                               "Sites_overlap_percentage",
                                               "Regions",
                                               "Regions_overlap_percentage"),
                                Coverage_thresholds, Diff_thresholds,
                                tsv_file = paste0(Path2, Peak_caller, "_iCLIPsites_overlap_with_Pos_Mutsites.tsv"))
      }
    }
  }
}


