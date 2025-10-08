## Set up working directory and library paths
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

## Load required libraries
library(Rsubread)
library(GenomicRanges)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(S4Vectors)
library(matrixStats)
library(dplyr)

################################################################################
## Function: group_sites_by_region
## Purpose: Aggregate mutation sites into contiguous genomic regions based on
##          proximity (group_thresh) and differential mutation significance
################################################################################
group_sites_by_region <- function(sites_paths, Coverage_thresholds, Diff_thresholds, group_thresh, filt_by) {
  for (i in seq_along(sites_paths)) {
    Res_path <- sites_paths[i]
    message("Processing path: ", Res_path)
    
    for (cov_thresh in Coverage_thresholds) {
      message("  Coverage threshold: ", cov_thresh)
      
      for (diff_thresh in Diff_thresholds) {
        message("    Diff threshold: ", diff_thresh)
        
        in_path <- paste0(Res_path, "cov_thresh_", cov_thresh, "/", filt_by, "/")
        sites_file <- paste0(in_path, "Pos_mut_thresh_", diff_thresh, ".tsv")
        region_file <- paste0(in_path, "Pos_mut_thresh_", diff_thresh, "_region.tsv")
        
        ## Read in mutation sites
        sites <- read.delim(sites_file, header = FALSE)
        colnames(sites) <- c("chr", "pos", "value")
        
        ## Group adjacent positions into regions
        grouped <- sites %>%
          select(chr, pos) %>%
          mutate(group = cumsum(diff(c(-Inf, pos)) >= group_thresh | chr != lag(chr, default = first(chr)))) %>%
          group_by(chr, group) %>%
          summarise(start = min(pos),
                    end = max(pos),
                    range = end - start + 1,
                    .groups = "drop") %>%
          select(chr, start, end, range)
        
        ## Write out grouped regions
        write.table(grouped, file = region_file, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
      }
    }
  }
}

################################################################################################
## Define mutation path structure
Mut_paths <- c(
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/"
)

ref_subfolder <- rep("ref_TotalRNA_cDNA_UV", length(Mut_paths))
subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")
group_thresh <- 5
Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
filt_by <- "FDR"

## Loop through all result paths and group mutation sites
for (sub1 in subfolders1) {
  for (sub2 in subfolders2) {
    sites_paths <- paste0(Mut_paths, sub1, "/", sub2, "/", ref_subfolder, "/")
    group_sites_by_region(sites_paths, Coverage_thresholds, Diff_thresholds, group_thresh, filt_by)
  }
}
