################################################################################
# mutation_stranded_annotation_pipeline.R
# Author: Mohamed B
# Description: Annotates mutation sites and regions by assigning strand info 
# from reference genome annotation. Outputs BED files for downstream analysis.
################################################################################

### Set working environment ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### Load libraries ###
library(Rsubread)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

################################################################################
# Function Definitions
################################################################################

# Load mutation sites and convert to GRanges with offset
load_sites <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = FALSE)
  Sites$strand <- "*"
  Sites <- Sites[!grepl("random", Sites$V1), ]
  GR <- GRanges(Sites$V1, IRanges(Sites$V2 - Left_offset, Sites$V2 + Right_offset), strand = Sites$strand)
  return(GR)
}

# Load mutation regions as GRanges
load_regions <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = TRUE)
  Sites$strand <- "*"
  Sites <- Sites[!is.na(Sites$start) & !grepl("random", Sites$chr), ]
  GR <- GRanges(Sites$chr, IRanges(Sites$start - Left_offset, Sites$end + Right_offset), strand = Sites$strand)
  return(GR)
}

# Load mutation region centers
load_regions_centre <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = TRUE)
  Sites$strand <- "*"
  Sites <- Sites[!is.na(Sites$start) & !grepl("random", Sites$chr), ]
  Sites$centre <- round((Sites$start + Sites$end)/2)
  GR <- GRanges(Sites$chr, IRanges(Sites$centre - Left_offset, Sites$centre + Right_offset), strand = Sites$strand)
  return(GR)
}

# Add strand info using genome annotation
add_strand <- function(Mut_sites_GR, genome_GR) {
  Overlaps <- findOverlaps(Mut_sites_GR, genome_GR)
  df <- data.frame(seqnames = seqnames(Mut_sites_GR)[queryHits(Overlaps)],
                   start = start(Mut_sites_GR)[queryHits(Overlaps)],
                   end = end(Mut_sites_GR)[queryHits(Overlaps)],
                   strand = strand(genome_GR)[subjectHits(Overlaps)])
  df <- unique(df)
  df$strand <- as.character(df$strand)
  
  # Collapse conflicting strand info
  df <- df %>%
    group_by(seqnames, start, end) %>%
    summarise(strand = if (n_distinct(strand) > 1) "*" else unique(strand), .groups = "drop") %>%
    filter(strand != "*")
  
  # BED formatting
  df$ID <- paste0("ID_", seq_len(nrow(df)))
  df$unused <- "*"
  df <- df[, c("seqnames", "start", "end", "ID", "unused", "strand")]
  
  return(df)
}

################################################################################
# Input Parameters
################################################################################

paths <- c("Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
           "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
           "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
           "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/")

ref_subfolder <- c("ref_P19_WT_TotalRNA_cDNA_no_UV",
                   "ref_P19_WT_TotalRNA_cDNA_no_UV",
                   "ref_TotalRNA_cDNA_UV",
                   "ref_TotalRNA_cDNA_UV")

subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
wins <- c(100)
filt_by <- "FDR"

gtf_file <- "genome/mm9/Annotation/gencode.vM1.annotation.gtf"

################################################################################
# Load Genome Annotation
################################################################################

genome_GR <- import(gtf_file)
genome_GR <- genome_GR[!grepl("_random", seqnames(genome_GR))]
genome_GR <- reduce(genome_GR)

################################################################################
# Main Annotation Loop
################################################################################

for (subfolder1 in subfolders1) {
  for (subfolder2 in subfolders2) {
    Mut_paths <- paste0(paths, subfolder1, "/", subfolder2, "/", ref_subfolder, "/")
    
    for (i in seq_along(Mut_paths)) {
      for (j in seq_along(Coverage_thresholds)) {
        Mut_path <- paste0(Mut_paths[i], "cov_thresh_", Coverage_thresholds[j], "/")
        print(Mut_path)
        
        for (k in seq_along(Diff_thresholds)) {
          print(Diff_thresholds[k])
          
          for (win in wins) {
            # Load GRanges for mutation types
            Pos_sites_GR <- load_sites(paste0(Mut_path, filt_by, "/", "Pos_mut_thresh_", Diff_thresholds[k], ".tsv"),
                                       win/2, win/2)
            
            Pos_regions_GR <- load_regions(paste0(Mut_path, filt_by, "/", "Pos_mut_thresh_", Diff_thresholds[k], "_region.tsv"),
                                           win/2, win/2)
            
            Pos_regions_centre_GR <- load_regions_centre(paste0(Mut_path, filt_by, "/", "Pos_mut_thresh_", Diff_thresholds[k], "_region.tsv"),
                                                         win/2, win/2)
            
            # Add strand info & filter
            df_sites <- add_strand(Pos_sites_GR, genome_GR)
            df_regions <- add_strand(Pos_regions_GR, genome_GR)
            df_regions_centre <- add_strand(Pos_regions_centre_GR, genome_GR)
            
            # Export as BED (gzipped)
            write.table(df_sites, gzfile(paste0(Mut_path, filt_by, "/", "Pos_mut_thresh_", Diff_thresholds[k], "_w_", win, ".bed.gz"), "w"),
                        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
            
            write.table(df_regions, gzfile(paste0(Mut_path, filt_by, "/", "Pos_mut_thresh_", Diff_thresholds[k], "_w_", win, "_region.bed.gz"), "w"),
                        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
            
            write.table(df_regions_centre, gzfile(paste0(Mut_path, filt_by, "/", "Pos_mut_thresh_", Diff_thresholds[k], "_w_", win, "_region_centre.bed.gz"), "w"),
                        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
          }
        }
      }
    }
  }
}
