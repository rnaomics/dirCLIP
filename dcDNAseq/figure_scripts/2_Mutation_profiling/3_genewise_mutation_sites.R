#!/usr/bin/env Rscript
# Gene-wise Significant Mutation Site Analysis
#
# Purpose:
# Summarize significant mutation sites and regions at the gene level
# and calculate overlap with iCLIP peaks.

############################################################
# Gene-wise mutation site overlap with iCLIP binding sites
############################################################

##############################
# Set working directory
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
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(dplyr)

############################################################
# Function: Load iCLIP peaks and convert to GRanges
############################################################
load_clip_tags <- function(iclip_path,
                           peak_caller,
                           left_offset,
                           right_offset) {
  
  ##########################################################
  # Load iCLIP BED file
  ##########################################################
  iclip_bed_df <- read.delim(
    iclip_path,
    header = FALSE
  )
  
  ##########################################################
  # Extend coordinates on positive strand
  ##########################################################
  positive_idx <- iclip_bed_df$V6 == "+"
  
  iclip_bed_df$V2[positive_idx] <-
    iclip_bed_df$V2[positive_idx] - left_offset
  
  iclip_bed_df$V3[positive_idx] <-
    iclip_bed_df$V3[positive_idx] + right_offset
  
  ##########################################################
  # Extend coordinates on negative strand
  ##########################################################
  negative_idx <- iclip_bed_df$V6 == "-"
  
  iclip_bed_df$V2[negative_idx] <-
    iclip_bed_df$V2[negative_idx] - right_offset
  
  iclip_bed_df$V3[negative_idx] <-
    iclip_bed_df$V3[negative_idx] + left_offset
  
  ##########################################################
  # Convert to GRanges
  ##########################################################
  iclip_gr <- GRanges(
    seqnames = iclip_bed_df[, "V1"],
    IRanges(
      iclip_bed_df[, "V2"],
      iclip_bed_df[, "V3"]
    ),
    iclip_bed_df[, "V6"]
  )
  
  return(iclip_gr)
}

############################################################
# Function: Load mutation sites and convert to GRanges
############################################################
load_mutation_sites <- function(path,
                                left_offset,
                                right_offset) {
  
  ##########################################################
  # Load mutation sites
  ##########################################################
  sites_df <- read.delim(
    path,
    header = FALSE
  )
  
  sites_df$strand <- "*"
  
  ##########################################################
  # Remove missing coordinates
  ##########################################################
  sites_df <- sites_df[
    !is.na(sites_df$V2),
  ]
  
  ##########################################################
  # Remove random chromosomes
  ##########################################################
  random_chr_idx <- grep(
    "random",
    sites_df$V1
  )
  
  if (!isEmpty(random_chr_idx)) {
    sites_df <- sites_df[-random_chr_idx, ]
  }
  
  ##########################################################
  # Convert to GRanges
  ##########################################################
  mutation_gr <- GRanges(
    seqnames = sites_df[, "V1"],
    IRanges(
      sites_df[, "V2"] - left_offset,
      sites_df[, "V2"] + right_offset
    ),
    sites_df$strand
  )
  
  return(mutation_gr)
}

############################################################
# Function: Load mutation regions and convert to GRanges
############################################################
load_mutation_regions <- function(path,
                                  left_offset,
                                  right_offset) {
  
  ##########################################################
  # Load mutation regions
  ##########################################################
  regions_df <- read.delim(
    path,
    header = TRUE
  )
  
  regions_df$strand <- "*"
  
  ##########################################################
  # Remove missing coordinates
  ##########################################################
  regions_df <- regions_df[
    !is.na(regions_df$start),
  ]
  
  ##########################################################
  # Remove random chromosomes
  ##########################################################
  random_chr_idx <- grep(
    "random",
    regions_df$chr
  )
  
  if (!isEmpty(random_chr_idx)) {
    regions_df <- regions_df[-random_chr_idx, ]
  }
  
  ##########################################################
  # Convert to GRanges
  ##########################################################
  mutation_gr <- GRanges(
    seqnames = regions_df[, "chr"],
    IRanges(
      regions_df[, "start"] - left_offset,
      regions_df[, "end"] + right_offset
    ),
    regions_df$strand
  )
  
  return(mutation_gr)
}

############################################################
# Function: Generate gene-level overlap summary
############################################################
summarize_gene_overlap <- function(mut_sites_gr,
                                   gene_annotation_gr,
                                   iclip_gr) {
  
  ##########################################################
  # Identify overlaps
  ##########################################################
  overlaps_with_gene <- findOverlaps(
    mut_sites_gr,
    gene_annotation_gr
  )
  
  overlaps_with_iclip <- findOverlaps(
    mut_sites_gr,
    iclip_gr
  )
  
  ##########################################################
  # Associate mutation sites with genes
  ##########################################################
  mutation_gene_df <- data.frame(
    mut_seqnames =
      seqnames(mut_sites_gr)[queryHits(overlaps_with_gene)],
    
    mut_position =
      start(mut_sites_gr)[queryHits(overlaps_with_gene)],
    
    gene_id =
      gene_annotation_gr$gene_id[
        subjectHits(overlaps_with_gene)
      ],
    
    Chrom =
      seqnames(gene_annotation_gr)[
        subjectHits(overlaps_with_gene)
      ],
    
    start =
      start(gene_annotation_gr)[
        subjectHits(overlaps_with_gene)
      ],
    
    end =
      end(gene_annotation_gr)[
        subjectHits(overlaps_with_gene)
      ]
  )
  
  ##########################################################
  # Mark iCLIP-overlapping mutation sites
  ##########################################################
  mutation_gene_df$Overlap_with_ICLIP <-
    ifelse(
      queryHits(overlaps_with_gene) %in%
        queryHits(overlaps_with_iclip),
      1,
      0
    )
  
  ##########################################################
  # Keep unique entries
  ##########################################################
  mutation_gene_df <- unique(mutation_gene_df)
  
  ##########################################################
  # Summarize overlap at gene level
  ##########################################################
  gene_summary_df <- mutation_gene_df %>%
    group_by(gene_id) %>%
    summarise(
      N_sites = n(),
      Sites_overlap_with_iCLIP =
        round(mean(Overlap_with_ICLIP), 2)
    )
  
  ##########################################################
  # Merge genomic coordinates and summary metrics
  ##########################################################
  gene_overlap_df <- mutation_gene_df %>%
    select(
      gene_id,
      Chrom,
      start,
      end
    ) %>%
    distinct() %>%
    left_join(
      gene_summary_df,
      by = "gene_id"
    )
  
  return(gene_overlap_df)
}

############################################################
# Function: Calculate gene-level overlap metrics for all
# peak callers
############################################################
get_genewise_mutation_overlap <- function(
    gene_annotation_gr,
    mut_sites_gr,
    peak_callers,
    iclip_paths) {
  
  for (peak_idx in seq_along(peak_callers)) {
    
    peak_caller <- peak_callers[peak_idx]
    iclip_path <- iclip_paths[peak_idx]
    
    iclip_gr <- load_clip_tags(
      iclip_path = iclip_path,
      peak_caller = peak_caller,
      left_offset = 0,
      right_offset = 0
    )
    
    overlap_df <- summarize_gene_overlap(
      mut_sites_gr = mut_sites_gr,
      gene_annotation_gr = gene_annotation_gr,
      iclip_gr = iclip_gr
    )
    
    if (peak_idx == 1) {
      
      gene_overlap_df <- overlap_df
      
    } else {
      
      gene_overlap_df <- cbind(
        gene_overlap_df,
        overlap_df[, 6]
      )
    }
  }
  
  colnames(gene_overlap_df)[-(1:5)] <-
    paste0(
      "Sites_overlap_with_",
      peak_callers
    )
  
  return(gene_overlap_df)
}

############################################################
# Function: Create output directory
############################################################
create_output_folder <- function(path) {
  
  ifelse(
    !dir.exists(path),
    dir.create(path, recursive = TRUE),
    FALSE
  )
}

############################################################
# Input mutation directories
############################################################
mutation_root_paths <- c(
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/",
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/"
)

reference_subfolders <- paste0(
  "ref_",
  c(
    "TotalRNA_cDNA_no_UV",
    "TotalRNA_cDNA_UV"
  )
)

analysis_subfolders <- c("2_Inc_del")

############################################################
# Peak caller settings
############################################################
peak_callers <- c("iCount")

iclip_paths <- c(
  "figure_scripts/figure_data/2_iclip/SRSF3/icount/SRSF3_merged.3nt_3nt.peaks.bed.gz"
)

############################################################
# Mutation filtering thresholds
############################################################
coverage_thresholds <- c(50)
difference_thresholds <- c(0.03)

############################################################
# Load gene annotations
############################################################
genome <- "mm9"
gtf <- "gencode"

gtf_file <- "genome/mm9/Annotation/gencode.vM1.annotation.gtf"

gene_annotation_gr <- import(gtf_file)

gene_annotation_gr <-
  gene_annotation_gr[
    gene_annotation_gr$type == "gene",
  ]

gene_annotation_gr$gene_id <-
  gene_annotation_gr$gene_name

############################################################
# Main analysis
############################################################
for (analysis_subfolder in analysis_subfolders) {
  
  mutation_input_paths <- paste0(
    mutation_root_paths,
    reference_subfolders,
    "/",
    analysis_subfolder,
    "/"
  )
  
  for (sample_idx in seq_along(mutation_input_paths)) {
    
    mutation_input_path <-
      mutation_input_paths[sample_idx]
    
    print(mutation_input_path)
    
    for (coverage_idx in seq_along(
      coverage_thresholds
    )) {
      
      coverage_threshold <-
        coverage_thresholds[coverage_idx]
      
      print(coverage_threshold)
      
      for (threshold_idx in seq_along(
        difference_thresholds
      )) {
        
        difference_threshold <-
          difference_thresholds[threshold_idx]
        
        print(difference_threshold)
        
        ####################################################
        # Output directories
        ####################################################
        threshold_output_path <- paste0(
          mutation_input_path,
          "cov_thresh_",
          coverage_threshold,
          "/"
        )
        
        annotation_output_path <- paste0(
          threshold_output_path,
          gtf,
          "/"
        )
        
        create_output_folder(
          annotation_output_path
        )
        
        ####################################################
        # Gene-wise mutation sites
        ####################################################
        mutation_gr <- load_mutation_sites(
          path = paste0(
            threshold_output_path,
            "Pos_mut_thresh_",
            difference_threshold,
            ".tsv"
          ),
          left_offset = 5,
          right_offset = 5
        )
        
        gene_overlap_df <-
          get_genewise_mutation_overlap(
            gene_annotation_gr =
              gene_annotation_gr,
            mut_sites_gr =
              mutation_gr,
            peak_callers =
              peak_callers,
            iclip_paths =
              iclip_paths
          )
        
        write.table(
          gene_overlap_df,
          paste0(
            annotation_output_path,
            "GW_mut_thresh_",
            difference_threshold,
            ".tsv"
          ),
          sep = "\t",
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
        
        ####################################################
        # Gene-wise mutation regions
        ####################################################
        mutation_gr <- load_mutation_regions(
          path = paste0(
            threshold_output_path,
            "Pos_mut_thresh_",
            difference_threshold,
            "_region.tsv"
          ),
          left_offset = 5,
          right_offset = 5
        )
        
        gene_overlap_df <-
          get_genewise_mutation_overlap(
            gene_annotation_gr =
              gene_annotation_gr,
            mut_sites_gr =
              mutation_gr,
            peak_callers =
              peak_callers,
            iclip_paths =
              iclip_paths
          )
        
        write.table(
          gene_overlap_df,
          paste0(
            annotation_output_path,
            "GW_mut_thresh_",
            difference_threshold,
            "_region.tsv"
          ),
          sep = "\t",
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
      }
    }
  }
}

############################################################
# End of script
############################################################