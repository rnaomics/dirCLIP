#!/usr/bin/env Rscript
# CNNC Motif Distance Analysis
#
# Purpose:
# Calculate positions of nearest CNNC motifs around mutation sites
# and compare against shuffled-sequence and GC-matched random controls.

# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

## Libraries ######
library(ShortRead)
library(stringr)
library(stringi)
library(ggplot2)
require(ggseqlogo) 
library(Biostrings)      
library(tidyverse)
library(writexl)
library(dplyr)

# ---------------------------------------------------------------------------------------------------
# 1. Get nearest CNNC from FASTA sequence
# ---------------------------------------------------------------------------------------------------
get_nearest_CNNC_from_fasta = function(
    fasta_path,
    Result_path,
    filename,
    kmers,
    CNNC_motifs,
    Length
){
  
  ####### Load FASTA ################################
  
  fasta_reads <- readFasta(
    dirPath = fasta_path
  )
  
  ####### Extract sequences #########################
  
  fasta_sequences =
    data.frame(
      fasta_reads@sread
    )[[1]]
  
  ####### Process each kmer #########################
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    ####### Find nearest CNNC motifs ################
    
    cnnc_locations =
      get_nearest_CNNC_locs(
        kmer,
        fasta_sequences,
        CNNC_motifs,
        Length
      )
    
    ####### Output ##################################
    
    output_path = paste0(
      Result_path,
      "kmer_",
      kmer,
      "/"
    )
    
    create_folder(output_path)
    
    write.table(
      cnnc_locations,
      paste0(output_path, filename),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  }
}

# ---------------------------------------------------------------------------------------------------
# 2. Get nearest CNNC from FASTA sequence after random shuffling
# ---------------------------------------------------------------------------------------------------

get_nearest_CNNC_from_fasta_after_shuffling = function(
    fasta_path,
    Result_path,
    randomer_filename,
    kmers,
    N,
    CNNC_motifs,
    Length
){
  
  ####### Load FASTA ################################
  
  fasta_reads <- readFasta(
    dirPath = fasta_path
  )
  
  ####### Extract sequences #########################
  
  fasta_sequences =
    data.frame(
      fasta_reads@sread
    )[[1]]
  
  ####### Process each kmer #########################
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    ####### Generate shuffled CNNC distances ########
    
    cnnc_locations =
      get_CNNC_motifs_dist_after_shuffling(
        fasta_sequences,
        N = N,
        kmer = kmer,
        CNNC_motifs,
        Length
      )
    
    ####### Output ##################################
    
    output_path = paste0(
      Result_path,
      "kmer_",
      kmer,
      "/"
    )
    
    create_folder(output_path)
    
    write.table(
      cnnc_locations,
      paste0(
        output_path,
        randomer_filename
      ),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  }
}

# ---------------------------------------------------------------------------------------------------
# 3. Get nearest CNNC from random background FASTA sequences
# ---------------------------------------------------------------------------------------------------

get_nearest_CNNC_from_fasta_rnd_bg_shuffling = function(
    fasta_path,
    Result_path,
    randomer_filename,
    kmers,
    N,
    CNNC_motifs,
    Length
){
  
  ####### Extract kmers #######
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    ####### Get CNNC locations #######
    
    cnnc_locations =
      get_CNNC_motifs_dist_rnd_bg_seq_shuffling(
        fasta_path,
        N = N,
        kmer = kmer,
        CNNC_motifs,
        Length
      )
    
    ####### Output #######
    
    output_path = paste0(
      Result_path,
      "kmer_",
      kmer,
      "/"
    )
    
    create_folder(output_path)
    
    write.table(
      cnnc_locations,
      paste0(
        output_path,
        randomer_filename
      ),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  }
}

# ---------------------------------------------------------------------------------------------------
# 4. Function to get CNNC motif distances from random background sequences
# ---------------------------------------------------------------------------------------------------

get_CNNC_motifs_dist_rnd_bg_seq_shuffling = function(
    fasta_path,
    N,
    kmer,
    CNNC_motifs,
    Length
){
  
  for (randomization_idx in 1:N){
    
    fasta_file = paste0(
      fasta_path,
      "rnd_bg_",
      randomization_idx,
      ".fasta.gz"
    )
    
    fasta_reads <- readFasta(
      dirPath = fasta_file
    )
    
    ####### Extract sequences ##########
    
    fasta_sequences =
      data.frame(
        fasta_reads@sread
      )[[1]]
    
    #### Randomly shuffle sequences ####
    fasta_sequences = stri_rand_shuffle(fasta_sequences)
    
    ####### Get CNNC locations #########
    
    cnnc_locations =
      get_nearest_CNNC_locs(
        kmer,
        fasta_sequences,
        CNNC_motifs,
        Length
      )
    
    if (randomization_idx == 1){
      
      cnnc_distance_df =
        cnnc_locations
      
    } else {
      
      cnnc_distance_df =
        merge(
          cnnc_distance_df,
          cnnc_locations,
          by = c("CNNC_locs"),
          all = TRUE
        )
    }
    
    print(randomization_idx)
  }
  
  ##### Replace NAs by 0 #####
  
  cnnc_distance_df[
    is.na(cnnc_distance_df)
  ] = 0
  
  ##### Rename frequency columns #####
  
  colnames(cnnc_distance_df)[-1] =
    paste0(
      "C",
      1:N,
      "_freq"
    )
  
  return(cnnc_distance_df)
}

# ---------------------------------------------------------------------------------------------------
# 5. Function to randomly shuffle sequences and calculate CNNC distances
# ---------------------------------------------------------------------------------------------------

get_CNNC_motifs_dist_after_shuffling = function(
    sequences,
    N,
    kmer,
    CNNC_motifs,
    Length
){
  
  for (randomization_idx in 1:N){
    
    #### Randomly shuffle sequences ####
    
    shuffled_sequences =
      stri_rand_shuffle(sequences)
    
    #### Get nearest CNNC locations ####
    
    cnnc_locations =
      get_nearest_CNNC_locs(
        kmer,
        shuffled_sequences,
        CNNC_motifs,
        Length
      )
    
    if (randomization_idx == 1){
      
      cnnc_distance_df =
        cnnc_locations
      
    } else {
      
      cnnc_distance_df =
        merge(
          cnnc_distance_df,
          cnnc_locations,
          by = c("CNNC_locs"),
          all = TRUE
        )
    }
    
    print(randomization_idx)
  }
  
  ##### Replace NAs by 0 #####
  
  cnnc_distance_df[
    is.na(cnnc_distance_df)
  ] = 0
  
  ##### Rename columns #####
  
  colnames(cnnc_distance_df)[-1] =
    paste0(
      "C",
      1:N,
      "_freq"
    )
  
  return(cnnc_distance_df)
}

# ---------------------------------------------------------------------------------------------------
# 6. Function to check CNNC
# ---------------------------------------------------------------------------------------------------

check_CNNC = function(motif) {
  
  k <- nchar(motif)
  
  any(
    sapply(
      1:(k - 3),
      function(i)
        substr(motif, i, i) == "C" &&
        substr(motif, i + 3, i + 3) == "C"
    )
  )
}

# ---------------------------------------------------------------------------------------------------
# 7. Function to get location of the nearest CNNC
# ---------------------------------------------------------------------------------------------------

get_nearest_CNNC = function(cnnc_matches){
  
  if (sum(cnnc_matches) >= 1){
    
    cnnc_positions =
      which(cnnc_matches) -
      length(cnnc_matches) / 2
    
    return(
      cnnc_positions[
        which.min(abs(cnnc_positions))
      ]
    )
    
  } else {
    
    return(Inf)
  }
}

# ---------------------------------------------------------------------------------------------------
# 8. Function to get nearest CNNC locations
# ---------------------------------------------------------------------------------------------------

get_nearest_CNNC_locs = function(
    kmer,
    sequences,
    CNNC_motifs,
    Length
){
  
  if (isEmpty(sequences)) {
    
    cnnc_locations =
      c(
        (1 - Length/2 + 2):
          (Length - Length/2 - 2),
        Inf
      )
    
    cnnc_locations =
      data.frame(
        CNNC_locs = cnnc_locations,
        Freq = 0
      )
  }
  
  else {
    
    sequence_lengths =
      nchar(sequences)
    
    unique_lengths =
      unique(sequence_lengths)
    
    n_possible_kmers =
      unique_lengths - kmer + 1
    
    motif_matrix =
      sapply(
        1:n_possible_kmers,
        function(x)
          substr(
            sequences,
            x,
            kmer + x - 1
          )
      )
    
    ##### Get CNNC motifs #####
    
    if (length(CNNC_motifs) == 0){
      
      unique_motifs =
        unique(
          as.vector(motif_matrix)
        )
      
      cnnc_motifs =
        sapply(
          unique_motifs,
          check_CNNC
        )
      
      cnnc_motifs =
        names(cnnc_motifs)[cnnc_motifs]
      
    } else {
      
      cnnc_motifs =
        CNNC_motifs
    }
    
    ##### Identify CNNC motifs #####
    
    cnnc_matches <- matrix(
      motif_matrix %in% cnnc_motifs,
      nrow = nrow(motif_matrix),
      ncol = ncol(motif_matrix)
    )
    
    ##### Get nearest CNNC location #####
    
    CNNC_locs = c()
    
    for (sequence_idx in 1:dim(cnnc_matches)[1]) {
      
      CNNC_locs[sequence_idx] =
        get_nearest_CNNC(
          cnnc_matches[sequence_idx, ]
        )
    }
    
    ##### Count frequencies #####
    
    full_distance_range =
      data.frame(
        CNNC_locs =
          c(
            (1 - Length/2 + 2):
              (Length - Length/2 - 2),
            Inf
          )
      )
    
    cnnc_frequency_df =
      data.frame(
        table(CNNC_locs)
      )

    cnnc_locations =
      merge(
        full_distance_range,
        cnnc_frequency_df,
        by = "CNNC_locs",
        all.x = TRUE
      )
    
    cnnc_locations =
      cnnc_locations %>%
      mutate(
        Freq = ifelse(
          is.na(Freq),
          0,
          Freq
        )
      )
  }
  
  return(cnnc_locations)
}

# ---------------------------------------------------------------------------------------------------
# 9. Function to Create folder
# ---------------------------------------------------------------------------------------------------

create_folder = function(Path){
  ifelse(!dir.exists(Path), dir.create(Path, recursive = TRUE), FALSE)
}

################################ SRSF3 Mg ##########################################################################################################

motif_analysis_paths = c(
  "Result/figure/2_Motif_analysis/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/"
)

# Top 10 iCLIP CNNC motifs (random background + shuffling)
CNNC_motifs = c(
  "TCTTC", "CTTCT", "CTTCC", "CTTCA", "CAGCA",
  "CTCCT", "CATCT", "TCCTC", "CCTTC", "CCTCC"
)

################## Setting thresholds ###########################################

Coverage_thresholds = c(50)
Diff_thresholds = c(0.03)

##################################################################################

kmers = c(5)

window_length = 100

mut_types = c(
  "sites"
)

##################################################################################

# CNNC_motifs = c() # Leave empty to identify all CNNC motifs

for (path_idx in seq_along(motif_analysis_paths)) {
  
  for (coverage_idx in seq_along(Coverage_thresholds)) {
    
    for (mut_type in mut_types) {
      
      motif_path = paste0(
        motif_analysis_paths[path_idx],
        "cov_thresh_",
        Coverage_thresholds[coverage_idx],
        "/",
        mut_type,
        "/win_", window_length, "/"
      )
      
      motif_rnd_bg_path = paste0(
        motif_analysis_paths[path_idx],
        "cov_thresh_",
        Coverage_thresholds[coverage_idx],
        "/",
        mut_type,
        "_rnd/win_", window_length, "/"
      )
      
      print(motif_path)
      
      for (threshold_idx in seq_along(Diff_thresholds)) {
        
        print(Diff_thresholds[threshold_idx])
        
        ########################################################
        # Input files
        ########################################################
        
        Pos_fasta_file = paste0(
          "Pos_mut_thresh_",
          Diff_thresholds[threshold_idx],
          ".fasta.gz"
        )
        
        CNNC_dist_file = paste0(
          "CNNC_dist_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        CNNC_dist_randomers_file = paste0(
          "CNNC_dist_randomers_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        CNNC_dist_rnd_bg_file = paste0(
          "CNNC_dist_rnd_bg_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        fasta_rnd_bg_region = paste0(
          "MT_",
          Diff_thresholds[threshold_idx],
          "/fasta/"
        )
        
        ########################################################
        # 1. Get nearest CNNC motifs from FASTA
        ########################################################
        
        print("Getting motifs from fasta...")
        
        get_nearest_CNNC_from_fasta(
          fasta_path = paste0(
            motif_path,
            Pos_fasta_file
          ),
          Result_path = motif_path,
          filename = CNNC_dist_file,
          kmers = kmers,
          CNNC_motifs,
          Length = window_length
        )
        
        ########################################################
        # 2. Generate shuffled controls
        ########################################################
        
        print("Getting random motifs from fasta...")
        
        N = 3
        set.seed(1)
        
        get_nearest_CNNC_from_fasta_after_shuffling(
          fasta_path = paste0(
            motif_path,
            Pos_fasta_file
          ),
          Result_path = motif_path,
          randomer_filename = CNNC_dist_randomers_file,
          kmers = kmers,
          N = N,
          CNNC_motifs,
          Length = window_length
        )
        
        ########################################################
        # 3. Generate random-background controls
        ########################################################
        
        print("Getting random motifs from fasta...")
        
        N = 3
        set.seed(1)
        
        get_nearest_CNNC_from_fasta_rnd_bg_shuffling(
          fasta_path = paste0(
            motif_rnd_bg_path,
            fasta_rnd_bg_region
          ),
          Result_path = motif_path,
          randomer_filename = CNNC_dist_rnd_bg_file,
          kmers = kmers,
          N = N,
          CNNC_motifs,
          Length = window_length
        )
      }
    }
  }
}
