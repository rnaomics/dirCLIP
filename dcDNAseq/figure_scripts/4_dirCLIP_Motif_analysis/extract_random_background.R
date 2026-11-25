#!/usr/bin/env Rscript

###############################################################
# GC-Matched Random Background Generation for Motif Analysis
#
# Author: Mohamed Bahrudeen
#
# Purpose:
# Generate GC- and region-matched random control sequences
# for downstream motif analysis.
#
# Inputs:
#   Pos_mut_thresh_*.tsv
#
# Outputs:
#   rnd_bg_*.bed.gz
#   rnd_bg_*.fasta.gz
#
###############################################################

# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(dplyr)
library(txdbmaker)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ggplot2)
library(rtracklayer)

genome <- BSgenome.Mmusculus.UCSC.mm9

# ------------------------------------------------------------
# 1. Function to load annotation
# ------------------------------------------------------------

load_annotation = function(annotation_gtf){
  
  gtf_granges <- import(annotation_gtf)
  
  genes_gr <- gtf_granges[
    gtf_granges$type == "gene"
  ]
  
  txdb_object <- makeTxDbFromGFF(annotation_gtf)
  
  exons <- exonsBy(
    txdb_object,
    by = "tx",
    use.names = TRUE
  )
  
  cds <- cdsBy(
    txdb_object,
    by = "tx",
    use.names = TRUE
  )
  
  utr5 <- fiveUTRsByTranscript(txdb_object)
  utr3 <- threeUTRsByTranscript(txdb_object)
  
  introns <- intronsByTranscript(
    txdb_object,
    use.names = TRUE
  )
  
  # Flatten features
  utr5 <- unlist(utr5)
  cds <- unlist(cds)
  utr3 <- unlist(utr3)
  exons <- unlist(exons)
  introns <- unlist(introns)
  
  annotation_list <- list(
    genes_gr = genes_gr,
    utr5 = utr5,
    cds = cds,
    utr3 = utr3,
    exons = exons,
    introns = introns
  )
  
  return(annotation_list)
}

# ---------------------------------------------------------------------------------------------------
# 2. Function to load dirCLIP sites, remove random chr, add strand information and annotate regions
# ---------------------------------------------------------------------------------------------------
annotate_region = function(
    dirCLIP_file,
    annotation_list
){
  
  genes_gr <- annotation_list$genes_gr
  utr5 <- annotation_list$utr5
  cds <- annotation_list$cds
  utr3 <- annotation_list$utr3
  exons <- annotation_list$exons
  introns <- annotation_list$introns
  
  ###### 1. Load dirCLIP sites ################
  
  dirclip_df <- read.table(
    dirCLIP_file,
    sep = "\t"
  )
  
  colnames(dirclip_df) <- c(
    "chr",
    "pos",
    "score"
  )
  
  ###### 2. Remove random chromosomes #########
  
  filtered_dirclip_df <- dirclip_df[
    !grepl("random", dirclip_df$chr),
  ]
  
  ###### 3. Create GRanges object #############
  
  site_granges <- GRanges(
    seqnames = filtered_dirclip_df$chr,
    ranges = IRanges(
      start = filtered_dirclip_df$pos,
      end = filtered_dirclip_df$pos
    )
  )
  
  ########## 4. Add strand information ########
  
  overlap_hits <- findOverlaps(
    site_granges,
    genes_gr
  )
  
  strand(site_granges)[queryHits(overlap_hits)] <-
    strand(genes_gr)[subjectHits(overlap_hits)]
  
  ########## 5. Annotate regions ##############
  
  region_annotation <- rep(
    "intergenic",
    length(site_granges)
  )
  
  region_annotation[
    countOverlaps(site_granges, introns) > 0
  ] <- "intron"
  
  region_annotation[
    countOverlaps(site_granges, exons) > 0
  ] <- "exon"
  
  region_annotation[
    countOverlaps(site_granges, cds) > 0
  ] <- "CDS"
  
  region_annotation[
    countOverlaps(site_granges, utr5) > 0
  ] <- "5UTR"
  
  region_annotation[
    countOverlaps(site_granges, utr3) > 0
  ] <- "3UTR"
  
  site_granges$region =
    region_annotation
  
  return(site_granges)
}

# ---------------------------------------------------------------------------------------------------
# 3. Function to generate GC-content pool of random regions 
# ---------------------------------------------------------------------------------------------------

make_gc_pool_random = function(gr, genome, win) {
  
  # 1. Sample uniform random positions inside each range
  
  random_fraction <- runif(length(gr))
  
  random_positions <-
    start(gr) +
    floor(
      (end(gr) - start(gr) + 1) *
        random_fraction
    )
  
  # 2. Create windows centered on sampled positions
  
  random_windows <- GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(
      start = random_positions - win %/% 2,
      width = win
    ),
    strand = strand(gr)
  )
  
  # 3. Get chromosome lengths
  
  chromosome_lengths <-
    seqlengths(genome)[
      as.character(
        seqnames(random_windows)
      )
    ]
  
  # 4. Remove invalid windows
  
  valid_window_idx <-
    start(random_windows) >= 1 &
    end(random_windows) <= chromosome_lengths
  
  random_windows <-
    random_windows[valid_window_idx]
  
  # 5. Extract sequence and compute GC content
  
  window_sequences <-
    getSeq(
      genome,
      random_windows
    )
  
  gc_content <-
    rowSums(
      letterFrequency(
        window_sequences,
        c("G", "C"),
        as.prob = TRUE
      )
    )
  
  mcols(random_windows)$gc <-
    gc_content
  
  random_windows
}

# ---------------------------------------------------------------------------------------------------
# 4. Function to get GC-matched random regions
# ---------------------------------------------------------------------------------------------------
get_GC_matched_random_region = function(
    gr,
    genome,
    GR_list,
    win = 100,
    max_gc_diff = 0.05,
    n_draws = 100
){
  
  start_time <- Sys.time()
  
  ###### 1. Precompute observed GC content #############################
  
  observed_windows <-
    resize(
      gr,
      width = win,
      fix = "center"
    )
  
  observed_sequences <-
    getSeq(
      genome,
      observed_windows
    )
  
  observed_gc <-
    calc_gc(observed_sequences)
  
  ###### 2. Generate random GC pools ##################################
  
  gc_matched_pools <- list(
    "5UTR" = make_gc_pool_random(
      GR_list$utr5,
      genome,
      win
    ),
    
    "CDS" = make_gc_pool_random(
      GR_list$cds,
      genome,
      win
    ),
    
    "3UTR" = make_gc_pool_random(
      GR_list$utr3,
      genome,
      win
    ),
    
    "exon" = make_gc_pool_random(
      GR_list$exons,
      genome,
      win
    ),
    
    "intron" = make_gc_pool_random(
      GR_list$introns,
      genome,
      win
    )
  )
  
  ########################################################
  
  random_windows_by_region <-
    vector(
      "list",
      length(unique(gr$region))
    )
  
  names(random_windows_by_region) <-
    unique(gr$region)
  
  for (region in unique(gr$region)) {
    
    region_idx <-
      which(gr$region == region)
    
    random_windows_by_region[[region]] <-
      process_region(
        gr_sub = gr[region_idx],
        real_gc = observed_gc[region_idx],
        pool = gc_matched_pools[[region]],
        region = region,
        max_gc_diff = max_gc_diff,
        n_draws = n_draws
      )
  }

  # random_windows <-
  #   do.call(
  #     c,
  #     unname(random_windows_by_region)
  #   )
  
  random_windows <- GRanges()
  for (region in names(random_windows_by_region)) {
    random_windows <- c(random_windows, random_windows_by_region[[region]])
  }
  
  end_time <- Sys.time()
  
  message(
    "GC-matched randomization runtime: ",
    round(
      difftime(
        end_time,
        start_time,
        units = "secs"
      ),
      2
    ),
    " seconds"
  )
  
  gc_difference <-
    observed_gc -
    mcols(random_windows)$gc
  
  mcols(random_windows)$gc_diff <-
    gc_difference
  
  return(random_windows)
}

# ---------------------------------------------------------------------------------------------------
# 5. Function to calculate GC
# ---------------------------------------------------------------------------------------------------
calc_gc <- function(sequences){
  
  rowSums(
    letterFrequency(
      sequences,
      letters = c("G", "C"),
      as.prob = TRUE
    )
  )
}

# ---------------------------------------------------------------------------------------------------
# 6. Function to extract windows
# ---------------------------------------------------------------------------------------------------

get_window <- function(
    gr,
    genome,
    window_width = 50
){
  
  window_size <- promoters(
    gr,
    upstream = window_width,
    downstream = window_width
  )
  
  return(window_size)
}

# ---------------------------------------------------------------------------------------------------
# 7. Region-wise processing
# ---------------------------------------------------------------------------------------------------

process_region <- function(
    gr_sub,
    real_gc,
    pool,
    region,
    max_gc_diff,
    n_draws = 20
){
  
  pool_gc_content <- mcols(pool)$gc
  
  n_observed <- length(real_gc)
  
  n_pool_regions <- length(pool_gc_content)
  
  sampled_indices <- matrix(
    sample.int(
      n_pool_regions,
      n_observed * n_draws,
      replace = TRUE
    ),
    nrow = n_observed
  )
  
  sampled_gc_content <- matrix(
    pool_gc_content[sampled_indices],
    nrow = n_observed
  )
  
  gc_difference <-
    abs(sampled_gc_content - real_gc)
  
  gc_match_idx <-
    gc_difference <= max_gc_diff
  
  # Position of minimum GC difference
  
  closest_gc_idx <-
    max.col(
      -gc_difference,
      ties.method = "first"
    )
  
  first_match_idx <-
    apply(
      gc_match_idx,
      1,
      function(x){
        
        idx <- which(x)
        
        if (length(idx) == 0){
          NA_integer_
        } else {
          idx[1]
        }
      }
    )
  
  selected_column <- ifelse(
    is.na(first_match_idx),
    closest_gc_idx,
    first_match_idx
  )
  
  selected_indices <-
    sampled_indices[
      cbind(
        seq_len(n_observed),
        selected_column
      )
    ]
  
  matched_regions <- pool[selected_indices]
  
  mcols(matched_regions)$region <- region
  
  matched_regions
}

# ---------------------------------------------------------------------------------------------------
# 8. Function to make region vs GC difference box plot
# ---------------------------------------------------------------------------------------------------

make_region_vs_difference_box_plot = function(random_regions_gr){
  
  gc_difference_df <- as.data.frame(
    mcols(random_regions_gr)
  )
  
  boxplot_figure =
    ggplot(
      gc_difference_df,
      aes(
        x = region,
        y = gc_diff
      )
    ) +
    
    geom_boxplot(
      outlier.size = 0.5
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    labs(
      x = "Genomic region",
      y = "GC difference (real − random)",
      title = "GC matching across genomic regions"
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  return(boxplot_figure)
}

# ---------------------------------------------------------------------------------------------------
# 9. Create folder
# ---------------------------------------------------------------------------------------------------

create_folder = function(folder_path){
  
  ifelse(
    !dir.exists(folder_path),
    dir.create(
      folder_path,
      recursive = TRUE
    ),
    FALSE
  )
}

################################# Setting dirCLIP paths ###############################

dirclip_paths <- c(
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/"
)

result_paths <- c(
  "Result/figure/2_Motif_analysis/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/"
)

cov_thresh <- 50
mut_thresh <- 0.03

win <- 100
n_draws <- 100

annotation_gtf <-
  "genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz"

########## 1. Load annotation #####################

annotation_list <-
  load_annotation(annotation_gtf)

#########################################################

for (path_idx in seq_along(dirclip_paths)) {
  
  dirCLIP_file <- paste0(
    dirclip_paths[path_idx],
    "cov_thresh_",
    cov_thresh,
    "/Pos_mut_thresh_",
    mut_thresh,
    ".tsv"
  )
  
  result_path <- result_paths[path_idx]
  
  ########## 2. Load dirCLIP sites and annotate regions ##########
  
  annotated_sites_gr <-
    annotate_region(
      dirCLIP_file,
      annotation_list
    )
  
  annotated_sites_gr <-
    annotated_sites_gr[
      !annotated_sites_gr$region == "intergenic"
    ]
  
  ######### 3. Create random controls ###########################
  
  random_seed <- 123
  set.seed(random_seed)
  
  n_random_controls <- 100   #  10
  
  for (control_idx in 1:n_random_controls) {
    
    ##### 1. GC‑matched random regions ########################
    
    random_regions_gr <-
      get_GC_matched_random_region(
        annotated_sites_gr,
        genome,
        annotation_list,
        win = win,
        max_gc_diff = 0.05,
        n_draws = n_draws
      )
    
    boxplot_figure <-
      make_region_vs_difference_box_plot(
        random_regions_gr
      )
    
    print(boxplot_figure)
    
    ##### 2. Output paths ####################################
    
    output_path <- paste0(
      result_path,
      "cov_thresh_", cov_thresh,
      "/sites_rnd/win_", win,
      "/MT_", mut_thresh, "/"
    )
    
    create_folder(
      paste0(output_path, "bed/")
    )
    
    create_folder(
      paste0(output_path, "fasta/")
    )
    
    ##### 3. Write BED #######################################
    
    random_bed_file <- paste0(
      output_path,
      "bed/rnd_bg_",
      control_idx,
      ".bed.gz"
    )
    
    export(
      random_regions_gr,
      con = random_bed_file,
      format = "bed"
    )
    
    ##### 4. Extract FASTA ###################################
    
    random_sequences <-
      getSeq(
        genome,
        random_regions_gr
      )
    
    names(random_sequences) <-
      paste0(
        "ID_",
        1:length(random_sequences)
      )
    
    random_fasta_file <- paste0(
      output_path,
      "fasta/rnd_bg_",
      control_idx,
      ".fasta.gz"
    )
    
    writeXStringSet(
      random_sequences,
      filepath = random_fasta_file,
      format = "fasta",
      compress = TRUE
    )
  }
}
