#!/usr/bin/env Rscript

# iCLIP Sequence Extraction and Random Background Generation
#
# Purpose:
# Extract iCLIP peak sequences and generate GC-matched random
# background regions for motif-enrichment analysis.
#

# Original workflow:
# 1. Load iCLIP peaks
# 2. Extend peaks to fixed-width windows
# 3. Annotate genomic regions
# 4. Remove intergenic peaks
# 5. Export BED and FASTA files
# 6. Generate GC-matched random regions
# 7. Export random BED and FASTA files

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

# ------------------------------------------------------------
# 2. Function to annotate genomic regions
# ------------------------------------------------------------

annotate_region = function(
    region_df,
    annotation_list
){
  
  genes_gr = annotation_list$genes_gr
  utr5 = annotation_list$utr5
  cds = annotation_list$cds
  utr3 = annotation_list$utr3
  exons = annotation_list$exons
  introns = annotation_list$introns
  
  ###### 1. Get centre position of region ####################
  
  region_df$center <-
    floor(
      (region_df$start + region_df$end) / 2
    )
  win = region_df$end[1] - region_df$start[1]
  ###### 2. Remove random chromosomes ########################
  
  filtered_region_df <- region_df[
    !grepl(
      "random",
      region_df$seqnames
    ),
  ]
  
  ###### 3. Create GRanges object ############################
  
  region_granges <- GRanges(
    seqnames = filtered_region_df$seqnames,
    ranges = IRanges(
      start = filtered_region_df$center,
      end = filtered_region_df$center
    ),
    strand = filtered_region_df$strand
  )
  
  ########## 4. Annotate genomic regions #####################
  
  region_annotation <- rep(
    "intergenic",
    length(region_granges)
  )
  
  region_annotation[
    countOverlaps(
      region_granges,
      introns
    ) > 0
  ] <- "intron"
  
  region_annotation[
    countOverlaps(
      region_granges,
      exons
    ) > 0
  ] <- "exon"
  
  region_annotation[
    countOverlaps(
      region_granges,
      cds
    ) > 0
  ] <- "CDS"
  
  region_annotation[
    countOverlaps(
      region_granges,
      utr5
    ) > 0
  ] <- "5UTR"
  
  region_annotation[
    countOverlaps(
      region_granges,
      utr3
    ) > 0
  ] <- "3UTR"
  
    region_granges$region =
    region_annotation
  
  
    region_granges <- resize(
      region_granges,
    width = win,
    fix = "center"
  )
  
  return(region_granges)
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
  
  random_windows <-
    do.call(
      c,
      unname(random_windows_by_region)
    )
  
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
# 8. Create folder
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

# ---------------------------------------------------------------------------------------------------
# 9. Function to get GC-matched random region BED/FASTA files
# ---------------------------------------------------------------------------------------------------

get_GC_matched_random_regions = function(
    gr,
    genome,
    GR_list,
    win,
    n_draws,
    Result_path,
    N
){
  
  for (random_region_idx in 1:N){
    
    ##### 1. Generate GC‑matched random regions ########################
    
    random_regions_gr =
      get_GC_matched_random_region(
        gr,
        genome,
        GR_list,
        win = win,
        max_gc_diff = 0.05,
        n_draws = n_draws
      )
    
    ##### 2. Output path ###############################################
    
    output_path = paste0(
      Result_path,
      "/random_region/"
    )
    
    create_folder(output_path)
    
    ##### 3. Write BED #################################################
    
    random_bed_file = paste0(
      output_path,
      "random_region_",
      random_region_idx,
      ".bed.gz"
    )
    
    export(
      random_regions_gr,
      con = random_bed_file,
      format = "bed"
    )
    
    ##### 4. Write FASTA ###############################################
    
    random_sequences <- getSeq(
      genome,
      random_regions_gr
    )
    
    names(random_sequences) <- paste0(
      "ID_",
      seq_along(random_sequences)
    )
    
    random_fasta_file = paste0(
      output_path,
      "random_region_",
      random_region_idx,
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

# ---------------------------------------------------------------------------------------------------
# 10. Function to load and extend iCLIP peaks
# ---------------------------------------------------------------------------------------------------
extend_iCLIP_regions <- function(
    iclip_path,
    Peak_caller,
    Left_offset,
    Right_offset,
    From_peak_centre = TRUE,
    genome,
    GR_list,
    Result_path,
    bed_file,
    fasta_file
) {
  
  ## ---------------- Load BED ----------------
  
  iclip_bed <- read.delim(
    iclip_path,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(iclip_bed) < 6) {
    stop("Input must be at least BED6 format")
  }
  
  ## ---------------- Remove random chromosomes ----------------
  
  valid_chr_idx <-
    !grepl("random", iclip_bed$V1)
  
  iclip_bed <-
    iclip_bed[valid_chr_idx, ]
  
  ## ---------------- Compute peak centre ----------------
  
  peak_center <-
    floor(
      (iclip_bed$V2 + iclip_bed$V3) / 2
    )
  
  ## ---------------- Extend regions ----------------
  
  plus_strand_idx  <- iclip_bed$V6 == "+"
  minus_strand_idx <- iclip_bed$V6 == "-"
  
  if (From_peak_centre) {
    
    iclip_bed$V2[plus_strand_idx] <-
      peak_center[plus_strand_idx] - Left_offset
    
    iclip_bed$V3[plus_strand_idx] <-
      peak_center[plus_strand_idx] + Right_offset
    
    iclip_bed$V2[minus_strand_idx] <-
      peak_center[minus_strand_idx] - Right_offset
    
    iclip_bed$V3[minus_strand_idx] <-
      peak_center[minus_strand_idx] + Left_offset
    
  } else {
    
    iclip_bed$V2[plus_strand_idx] <-
      iclip_bed$V2[plus_strand_idx] - Left_offset
    
    iclip_bed$V3[plus_strand_idx] <-
      iclip_bed$V3[plus_strand_idx] + Right_offset
    
    iclip_bed$V2[minus_strand_idx] <-
      iclip_bed$V2[minus_strand_idx] - Right_offset
    
    iclip_bed$V3[minus_strand_idx] <-
      iclip_bed$V3[minus_strand_idx] + Left_offset
  }
  
  ## ---------------- Fix invalid coordinates ----------------
  
  chr_len <-
    seqlengths(genome)[iclip_bed$V1]
  
  valid_region_idx <-
    iclip_bed$V2 >= 0 &
    iclip_bed$V3 > iclip_bed$V2 &
    iclip_bed$V3 <= chr_len
  
  iclip_bed <-
    iclip_bed[valid_region_idx, ]
  
  ## ---------------- Rebuild BED6 ----------------
  
  iclip_bed$V4 <-
    paste0(
      "ID_",
      seq_len(nrow(iclip_bed))
    )
  
  iclip_bed$V5 <- "*"
  
  iclip_bed <-
    iclip_bed[, 1:6]
  
  colnames(iclip_bed) <-
    c(
      "seqnames",
      "start",
      "end",
      "ID",
      "unused",
      "strand"
    )
  
  ## ---------------- Annotate regions ----------------
  
  annotated_peaks_gr <-
    annotate_region(
      iclip_bed,
      GR_list
    )
  
  non_intergenic_idx <-
    annotated_peaks_gr$region != "intergenic"
  
  annotated_peaks_gr <-
    annotated_peaks_gr[
      non_intergenic_idx
    ]
  
  iclip_bed <-
    iclip_bed[
      non_intergenic_idx,
    ]
  
  ## ---------------- Write BED.gz ----------------
  
  dir.create(
    Result_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  con <- gzfile(
    file.path(Result_path, bed_file),
    "w"
  )
  
  write.table(
    iclip_bed,
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  close(con)
  
  ## ---------------- Write FASTA.gz ----------------
  
  peak_sequences <- getSeq(
    genome,
    annotated_peaks_gr
  )
  
  names(peak_sequences) <-
    paste0(
      "ID_",
      seq_along(peak_sequences)
    )
  
  fasta_output_file <-
    paste0(
      Result_path,
      fasta_file
    )
  
  writeXStringSet(
    peak_sequences,
    filepath = fasta_output_file,
    format = "fasta",
    compress = TRUE
  )

  return(annotated_peaks_gr)
}

############################### From peak centre ############################################

iclip_peak_paths = c(
  "figure_scripts/figure_data/2_iclip/SRSF3/icount/SRSF3_merged.3nt_3nt.peaks.bed.gz"
)

Peak_callers = c("icount")

result_paths = c(
  "Result/figure/7_iclip/SRSF3/icount/Motif_analysis/"
)

win = 100

From_peak_centre = TRUE

##################### Generate background BED file ##########################################

Generate_background_bed = TRUE

n_random_controls = 10   # 100

options(scipen = 999)

########### GENCODE annotation ##############################################################

gtf_file =
  "genome/mm9/Annotation/gencode.vM1.annotation.gtf"

########## 1. Load annotation #####################

annotation_list =
  load_annotation(gtf_file)

################################################################################################

for (peak_caller_idx in seq_along(Peak_callers)) {
  
  iclip_peak_path =
    iclip_peak_paths[peak_caller_idx]
  
  Result_path = paste0(
    result_paths[peak_caller_idx],
    "win_",
    win,
    "/"
  )
  
  ################ 1. Extend peaks and export BED/FASTA ####################################
  
  iclip_regions_gr =
    extend_iCLIP_regions(
      iclip_path = iclip_peak_path,
      Peak_caller = Peak_callers[peak_caller_idx],
      Left_offset = win / 2,
      Right_offset = win / 2,
      From_peak_centre = From_peak_centre,
      genome = genome,
      GR_list = annotation_list,
      Result_path = Result_path,
      bed_file = "iCLIP.bed.gz",
      fasta_file = "iCLIP.fasta.gz"
    )
  
  ################ 2. Generate GC-matched random controls ###################################
  
  set.seed(1)
  
  if (Generate_background_bed) {
    
    n_draws = 100
    
    get_GC_matched_random_regions(
      gr = iclip_regions_gr,
      genome = genome,
      GR_list = annotation_list,
      win = win,
      n_draws = n_draws,
      Result_path = Result_path,
      N = n_random_controls
    )
  }
}