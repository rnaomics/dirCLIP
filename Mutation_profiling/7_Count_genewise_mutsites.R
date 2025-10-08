################################################################################
# genewise_mutation_overlap_analysis.R
# Author: Mohamed B
# Description:
# This script reads Nanopore-derived mutation site and region data, compares
# genewise overlaps with iCLIP peaks from multiple peak callers, and outputs
# annotated tables for downstream analysis.
################################################################################

### Setup ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### Load Required Libraries ###
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
library(rtracklayer)
library(dplyr)

### Load External Functions ###
source("Tools/DIRCLIP_analysis/8_Mutation_profiling/1_mm9_genome/count_genewise_sites_func.R")

################################################################################
# Functions
################################################################################

load_CLIP_tags <- function(iclip_path, Peak_caller, Left_offset, Right_offset) {
  if (Peak_caller == "Michael_unknown") {
    iclip_bed <- read.table(iclip_path, skip = 1)
    iclip_bed$Strand <- ifelse(iclip_bed$V4 < 0, "-", "+")
    iclip_bed$V2[iclip_bed$Strand == "+"] <- iclip_bed$V2[iclip_bed$Strand == "+"] - Left_offset
    iclip_bed$V3[iclip_bed$Strand == "+"] <- iclip_bed$V3[iclip_bed$Strand == "+"] + Right_offset
    iclip_bed$V2[iclip_bed$Strand == "-"] <- iclip_bed$V2[iclip_bed$Strand == "-"] - Right_offset
    iclip_bed$V3[iclip_bed$Strand == "-"] <- iclip_bed$V3[iclip_bed$Strand == "-"] + Left_offset
    return(GRanges(iclip_bed$V1, IRanges(iclip_bed$V2, iclip_bed$V3), strand = iclip_bed$Strand))
  } else {
    iclip_bed <- read.delim(iclip_path, header = FALSE)
    iclip_bed$V2[iclip_bed$V6 == "+"] <- iclip_bed$V2[iclip_bed$V6 == "+"] - Left_offset
    iclip_bed$V3[iclip_bed$V6 == "+"] <- iclip_bed$V3[iclip_bed$V6 == "+"] + Right_offset
    iclip_bed$V2[iclip_bed$V6 == "-"] <- iclip_bed$V2[iclip_bed$V6 == "-"] - Right_offset
    iclip_bed$V3[iclip_bed$V6 == "-"] <- iclip_bed$V3[iclip_bed$V6 == "-"] + Left_offset
    return(GRanges(iclip_bed$V1, IRanges(iclip_bed$V2, iclip_bed$V3), strand = iclip_bed$V6))
  }
}

load_sites <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = FALSE)
  Sites$strand <- "*"
  Sites <- Sites[!grepl("random", Sites$V1), ]
  return(GRanges(Sites$V1, IRanges(Sites$V2 - Left_offset, Sites$V2 + Right_offset), strand = Sites$strand))
}

load_regions <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = TRUE)
  Sites$strand <- "*"
  Sites <- Sites[!is.na(Sites$start) & !grepl("random", Sites$chr), ]
  return(GRanges(Sites$chr, IRanges(Sites$start - Left_offset, Sites$end + Right_offset), strand = Sites$strand))
}

get_genewise_mut_sites4 <- function(gtf_GR, mut_sites_GR, Peak_callers, iclip_paths) {
  for (z in seq_along(Peak_callers)) {
    icLIP_GR <- load_CLIP_tags(iclip_paths[z], Peak_callers[z], 0, 0)
    df <- get_genewise_mut_sites3(gtf_GR = gtf_GR, mut_sites_GR = mut_sites_GR, icLIP_GR = icLIP_GR)
    df_gw <- if (z == 1) df else cbind(df_gw, df[,6])
  }
  colnames(df_gw)[-(1:5)] <- paste0("Sites_overlap_with_", Peak_callers)
  return(df_gw)
}

################################################################################
# Input Configuration
################################################################################

Mut_paths <- c("Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
               "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/")
ref_subfolder <- c("ref_TotalRNA_cDNA_UV", "ref_TotalRNA_cDNA_UV")
subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

Peak_callers <- c("CLIPper", "iCount", "pureclip", "piranha", "paraclu")
iclip_paths <- c("Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/clipper/SRSF3_merged.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/icount/SRSF3_merged.3nt_3nt.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/pureclip/SRSF3_merged.8nt.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/piranha/SRSF3_merged.3nt_3nt.peaks.bed.gz",
                 "Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/paraclu/SRSF3_merged.10_200nt_2.peaks.bed.gz")

gtf_file <- "genome/mm9/Annotation/gencode.vM1.annotation.gtf"
gtf_GR <- import(gtf_file)
gtf_GR <- gtf_GR[gtf_GR$type == "gene", ]
gtf_GR$gene_id <- gtf_GR$gene_name

Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
filt_by <- "FDR"

#################################################################################################################################
# Iterates over mutation datasets, thresholds, and conditions to generate
# genewise summaries of mutation site and region overlaps with iCLIP peaks.
################################################################################

for (subfolder1 in subfolders1) {
  for (subfolder2 in subfolders2) {
    
    # Construct mutation site paths for all replicates and conditions
    Mut_sites_paths <- paste0(Mut_paths, subfolder1, "/", subfolder2, "/", ref_subfolder, "/")
    
    for (i in seq_along(Mut_sites_paths)) {
      Mut_sites_path <- Mut_sites_paths[i]
      Res_path <- Mut_sites_path  # Same as input path; reuse
      
      cat("🔍 Processing:", Mut_sites_path, "\n")
      
      for (cov_thresh in Coverage_thresholds) {
        cat("  ➤ Coverage threshold:", cov_thresh, "\n")
        
        for (diff_thresh in Diff_thresholds) {
          cat("    ➤ Difference threshold:", diff_thresh, "\n")
          
          # Define result directory containing filtered mutation data
          res_path <- paste0(Res_path, "cov_thresh_", cov_thresh, "/", filt_by, "/")
          
          #### (1) Genewise mutation site overlap analysis ####
          
          # Load mutation sites as GRanges object
          mut_sites_GR <- load_sites(
            Path = paste0(res_path, "Pos_mut_thresh_", diff_thresh, ".tsv"),
            Left_offset = 5, Right_offset = 5
          )
          
          # Compute overlaps across all iCLIP peak callers
          df_gw_mut_sites <- get_genewise_mut_sites4(
            gtf_GR = gtf_GR,
            mut_sites_GR = mut_sites_GR,
            Peak_callers = Peak_callers,
            iclip_paths = iclip_paths
          )
          
          # Save genewise mutation site summary
          write.table(
            df_gw_mut_sites,
            file = paste0(res_path, "GW_mut_thresh_", diff_thresh, ".tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
          )
          
          #### (2) Genewise mutation region overlap analysis ####
          
          # Load mutation regions as GRanges
          mut_sites_GR <- load_regions(
            Path = paste0(res_path, "Pos_mut_thresh_", diff_thresh, "_region.tsv"),
            Left_offset = 5, Right_offset = 5
          )
          
          # Compute overlaps for regions
          df_gw_mut_sites <- get_genewise_mut_sites4(
            gtf_GR = gtf_GR,
            mut_sites_GR = mut_sites_GR,
            Peak_callers = Peak_callers,
            iclip_paths = iclip_paths
          )
          
          # Save genewise region summary
          write.table(
            df_gw_mut_sites,
            file = paste0(res_path, "GW_mut_thresh_", diff_thresh, "_region.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
          )
        }
      }
    }
  }
}
########################################################################################