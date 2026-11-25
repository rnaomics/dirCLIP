#!/usr/bin/env Rscript

# Motif Enrichment Analysis
# Author: Mohamed Bahrudeen
#
# Purpose:
# Identify enriched sequence motifs from significant mutation sites.
#
# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

## Libraries ######
library(ShortRead)
library(stringr)
library(stringi)
library(ggplot2)
# library(motifStack)
require(ggseqlogo) 
library(Biostrings)
library(tidyverse)
library(writexl)
library(dplyr)
library(matrixStats)

######################## Create folder ######################################################

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

####################### Function to get motifs from random background sequences ####################################

get_motifs_from_rnd_bg_seq = function(
    rnd_fasta_path,
    Result_path,
    motif_filename,
    randomer_filename,
    kmers,
    N
){
  
  ####### Extracting kmers #######
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    ####### Load discovered motifs #######
    
    top_motifs = read.delim(
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/",
        motif_filename
      )
    )
    
    ####### Get motif counts from random sequences #######
    
    random_motif_counts =
      get_rnd_bg_seq_motif_counts(
        rnd_fasta_path,
        Top_motifs = top_motifs,
        N = N,
        kmer = kmer
      )
    
    ####### Output #######
    
    output_path =
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/"
      )
    
    create_folder(output_path)
    
    write.table(
      random_motif_counts,
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

############################ Function to get motif counts from random background sequences #######################

get_rnd_bg_seq_motif_counts = function(
    rnd_fasta_path,
    Top_motifs,
    N,
    kmer
){
  
  top_motifs = Top_motifs[, -3]
  
  motif_reference = top_motifs
  
  for (random_idx in 1:N) {
    
    fasta_path = paste0(
      rnd_fasta_path,
      "rnd_bg_",
      random_idx,
      ".fasta.gz"
    )
    
    ####### Load FASTA #######
    
    fasta_reads <- readFasta(
      dirPath = fasta_path
    )
    
    ####### Extract sequences #######
    
    random_sequences =
      data.frame(
        fasta_reads@sread
      )[[1]]
    
    ####### Count motifs #######
    
    random_motif_counts =
      find_motifs(
        kmer,
        random_sequences
      )
    
    ####### Add missing motifs #######
    
    missing_motifs <- setdiff(
      motif_reference$motifs,
      random_motif_counts$motifs
    )
    
    if (length(missing_motifs) > 0) {
      
      random_motif_counts <- rbind(
        random_motif_counts,
        data.frame(
          motifs = missing_motifs,
          total_count = 0,
          prop = 0
        )
      )
    }
    
    random_motif_counts =
      random_motif_counts[, -3]
    
    ####### Merge counts #######
    
    top_motifs =
      merge(
        top_motifs,
        random_motif_counts,
        by.x = "motifs",
        by.y = "motifs"
      )
    
    print(random_idx)
  }
  
  ####### Rename columns #######
  
  column_names =
    c(
      "T",
      paste0("C", 1:N)
    )
  
  colnames(top_motifs)[-1] =
    column_names
  
  ####### Sort by observed motif counts #######
  
  top_motifs =
    top_motifs %>%
    arrange(desc(T))
  
  ####### Remove observed count column #######
  
  top_motifs =
    top_motifs[, -2]
  
  return(top_motifs)
}

####################### Get motifs after random shuffling ####################################

get_motifs_after_shuffling = function(
    fasta_path,
    Result_path,
    motif_filename,
    randomer_filename,
    kmers,
    N
){
  
  ## N = Number of random shufflings
  
  ####### Loading the FASTA file #########
  
  fasta_reads <- readFasta(
    dirPath = fasta_path
  )
  
  ####### Extract sequences ##########
  
  fasta_sequences =
    data.frame(
      fasta_reads@sread
    )[[1]]
  
  ####### Extract kmers #######
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    ####### Load motifs #######
    
    top_motifs = read.delim(
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/",
        motif_filename
      )
    )
    
    ####### Get motif counts after shuffling #######
    
    random_motif_counts =
      get_random_motif_counts(
        fasta_sequences,
        Top_motifs = top_motifs,
        N = N,
        kmer = kmer
      )
    
    ####### Write results #######
    
    output_path =
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/"
      )
    
    create_folder(output_path)
    
    write.table(
      random_motif_counts,
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

####################### Function to get enriched motifs by shuffled sequences as control ############################

get_enriched_motifs_by_shuffling = function(
    motif_path,
    Pos_fasta_file,
    Pos_motifs_file,
    Pos_randomers_file,
    Pos_zscore_file,
    Pos_seqlogo_file,
    kmers,
    N,
    z_thresh,
    pval_thresh
){
  
  ################ Get the motifs after random shuffling from the read sequences #################################
  
  print("Getting random motifs from fasta...")
  
  # N = 100 # number of runs
  set.seed(1)
  
  get_motifs_after_shuffling(
    fasta_path = paste0(
      motif_path,
      Pos_fasta_file
    ),
    Result_path = motif_path,
    motif_filename = Pos_motifs_file,
    randomer_filename = Pos_randomers_file,
    kmers = kmers,
    N = N
  )
  
  ###################### 3. Calculate the z-score ##################################################################
  
  print("Calculate the z-score...")
  
  get_motifs_z_score(
    Result_path = motif_path,
    motif_filename = Pos_motifs_file,
    randomer_filename = Pos_randomers_file,
    zscore_filename = Pos_zscore_file,
    kmers = kmers
  )
  
  ###################### 4. Calculate right-tailed Fisher's exact test ##############################################
  
  print("Calculate the right-tailed Fisher's exact test...")
  
  get_motifs_pvalue(
    Result_path = motif_path,
    zscore_filename = Pos_zscore_file,
    kmers = kmers
  )
  
  ##################### 5. Generate the sequence logo plot ##########################################################
  
  print("Get the sequence logo plot...")
  
  make_motifs_sequence_logo_plot(
    Result_path = motif_path,
    kmers = kmers,
    zscore_filename = Pos_zscore_file,
    seq_logo_filename = Pos_seqlogo_file,
    z_thresh = z_thresh,
    pval_thresh = pval_thresh
  )
}

####################### Function to get enriched motifs using random background sequences as control ############################

get_enriched_motifs_using_rnd_bg = function(
    motif_path,
    Pos_motifs_file,
    Pos_rnd_bg_file,
    Pos_zscore_file,
    Pos_seqlogo_file,
    motif_rnd_bg_path,
    fasta_rnd_bg_region,
    kmers,
    N,
    z_thresh,
    pval_thresh
){
  
  ################ Get motifs from random background sequences #################################
  
  print("Getting random motifs from fasta...")
  
  set.seed(1)
  
  get_motifs_from_rnd_bg_seq(
    rnd_fasta_path = paste0(
      motif_rnd_bg_path,
      fasta_rnd_bg_region
    ),
    Result_path = motif_path,
    motif_filename = Pos_motifs_file,
    randomer_filename = Pos_rnd_bg_file,
    kmers = kmers,
    N = N
  )
  
  ###################### Calculate z-scores ###################################################
  
  print("Calculate the z-score...")
  
  get_motifs_z_score(
    Result_path = motif_path,
    motif_filename = Pos_motifs_file,
    randomer_filename = Pos_rnd_bg_file,
    zscore_filename = Pos_zscore_file,
    kmers = kmers
  )
  
  ###################### Calculate right-tailed Fisher's exact test ############################
  
  print("Calculate the right-tailed Fisher's exact test...")
  
  get_motifs_pvalue(
    Result_path = motif_path,
    zscore_filename = Pos_zscore_file,
    kmers = kmers
  )
  
  ###################### Generate sequence logo plots ##########################################
  
  print("Get the sequence logo plot...")
  
  make_motifs_sequence_logo_plot(
    Result_path = motif_path,
    kmers = kmers,
    zscore_filename = Pos_zscore_file,
    seq_logo_filename = Pos_seqlogo_file,
    z_thresh = z_thresh,
    pval_thresh = pval_thresh
  )
}

############## Get motifs from FASTA sequences ###################################################

get_motifs_from_fasta = function(
    fasta_path,
    Result_path,
    filename,
    kmers
){
  
  ####### Load FASTA ############################
  
  fasta_reads <- readFasta(
    dirPath = fasta_path
  )
  
  ####### Extract sequences #####################
  
  fasta_sequences =
    data.frame(
      fasta_reads@sread
    )[[1]]
  
  ####### Extract kmers #########################
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    ####### Get motif counts ####################
    
    motif_counts =
      find_motifs(
        kmer,
        fasta_sequences
      )
    
    ####### Output ##############################
    
    output_path =
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/"
      )
    
    create_folder(output_path)
    
    write.table(
      motif_counts,
      paste0(output_path, filename),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  }
}

############################ Function to find the most frequent motifs ##########################

find_motifs = function(
    kmer,
    sequences
){
  
  if (isEmpty(sequences)) {
    
    motif_counts_df =
      data.frame(
        motifs = "Empty",
        total_count = 0,
        prop = 0
      )
    
  } else {
    
    sequence_lengths =
      nchar(sequences)
    
    unique_lengths =
      unique(sequence_lengths)
    
    all_kmers = c()
    
    for (sequence_length in unique_lengths) {
      
      length_matched_sequences =
        sequences[
          sequence_lengths == sequence_length
        ]
      
      n_possible_kmers =
        sequence_length - kmer + 1
      
      kmer_matrix =
        sapply(
          1:n_possible_kmers,
          function(x)
            substr(
              length_matched_sequences,
              x,
              kmer + x - 1
            )
        )
      
      all_kmers =
        c(
          all_kmers,
          as.vector(kmer_matrix)
        )
    }
    
    motif_counts_df =
      data.frame(
        motifs = all_kmers
      )
    
    motif_counts_df =
      motif_counts_df %>%
      group_by(motifs) %>%
      summarise(
        total_count = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(total_count))
    
    motif_counts_df$prop =
      motif_counts_df$total_count /
      sum(motif_counts_df$total_count)
  }
  
  return(motif_counts_df)
}

############################ Function to randomly shuffle sequences and count motifs #######################

get_random_motif_counts = function(
    sequences,
    Top_motifs,
    N,
    kmer
){
  
  top_motifs = Top_motifs[, -3]
  
  reference_motifs = top_motifs
  
  for (randomization_idx in 1:N) {
    
    #### Shuffle sequences #####################
    
    shuffled_sequences =
      stri_rand_shuffle(sequences)
    
    #### Count motifs ##########################
    
    random_motif_counts =
      find_motifs(
        kmer,
        shuffled_sequences
      )
    
    #### Add missing motifs ####################
    
    missing_motifs <- setdiff(
      reference_motifs$motifs,
      random_motif_counts$motifs
    )
    
    if (length(missing_motifs) > 0) {
      
      random_motif_counts <- rbind(
        random_motif_counts,
        data.frame(
          motifs = missing_motifs,
          total_count = 0,
          prop = 0
        )
      )
    }
    
    random_motif_counts =
      random_motif_counts[, -3]
    
    #### Merge counts ##########################
    
    top_motifs =
      merge(
        top_motifs,
        random_motif_counts,
        by.x = "motifs",
        by.y = "motifs"
      )
    
    print(randomization_idx)
  }
  
  #### Rename columns ##########################
  
  column_names =
    c(
      "T",
      paste0("C", 1:N)
    )
  
  colnames(top_motifs)[-1] =
    column_names
  
  #### Sort by observed frequency ##############
  
  top_motifs =
    top_motifs %>%
    arrange(desc(T))
  
  #### Remove observed-count column ###########
  
  top_motifs =
    top_motifs[, -2]
  
  return(top_motifs)
}

###################### Get motif z-scores ##################################################

get_motifs_z_score = function(
    Result_path,
    motif_filename,
    randomer_filename,
    zscore_filename,
    kmers
){
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    output_path =
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/"
      )
    
    ######## Load results ######################
    
    motif_counts <- read.delim(
      paste0(output_path, motif_filename)
    )
    
    random_motif_counts <- read.delim(
      paste0(output_path, randomer_filename)
    )
    
    ######## Calculate random means and SDs ####
    
    random_mean_counts <-
      rowMeans(
        random_motif_counts[, -1]
      )
    
    random_sd_counts <-
      rowSds(
        data.matrix(
          random_motif_counts[, -1]
        )
      )
    
    ######## Calculate z-scores ################
    
    z_scores =
      (
        motif_counts$total_count -
          random_mean_counts
      ) /
      random_sd_counts
    
    ######## Add statistics ###################
    
    motif_counts$shuffled_mean =
      random_mean_counts
    
    motif_counts$shuffled_std =
      random_sd_counts
    
    motif_counts$z_score =
      z_scores
    
    ######## Sort by z-score ##################
    
    motif_counts =
      motif_counts %>%
      arrange(desc(z_score))
    
    ######## Output ###########################
    
    create_folder(output_path)
    
    write.table(
      motif_counts,
      paste0(
        output_path,
        zscore_filename
      ),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  }
}

###################### Get Fisher's exact test p-values #####################################

get_motifs_pvalue = function(
    Result_path,
    zscore_filename,
    kmers
){
  
  ####################################################
  # Perform right-tailed Fisher's exact test
  ####################################################
  
  perform_fisher_test <- function(row) {
    
    contingency_table <- matrix(
      row,
      nrow = 2,
      byrow = TRUE
    )
    
    fisher_result <- fisher.test(
      contingency_table,
      alternative = "greater"
    )
    
    return(fisher_result$p.value)
  }
  
  ####################################################
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    output_path =
      paste0(
        Result_path,
        "kmer_",
        kmer,
        "/"
      )
    
    ######## Load z-score table ################
    
    motif_counts <- read.delim(
      paste0(
        output_path,
        zscore_filename
      )
    )
    
    ######## Build Fisher test table ##########
    
    total_observed_counts =
      sum(motif_counts$total_count)
    
    total_shuffled_counts =
      sum(motif_counts$shuffled_mean)
    
    fisher_input_df =
      motif_counts[
        ,
        c(
          "total_count",
          "shuffled_mean"
        )
      ]
    
    fisher_input_df$sub_tot_count =
      total_observed_counts -
      fisher_input_df$total_count
    
    fisher_input_df$sub_shuffle =
      total_shuffled_counts -
      fisher_input_df$shuffled_mean
    
    ######## Calculate p-values ###############
    
    motif_counts$p_value =
      apply(
        fisher_input_df,
        1,
        perform_fisher_test
      )
    
    ######## Sort by significance #############
    
    motif_counts =
      motif_counts %>%
      arrange(p_value)
    
    ######## Output ###########################
    
    create_folder(output_path)
    
    write.table(
      motif_counts,
      paste0(
        output_path,
        zscore_filename
      ),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
  }
}


####################### Making sequence logo plots using motifs above the z-score threshold ######################

make_motifs_sequence_logo_plot = function(
    Result_path,
    zscore_filename,
    seq_logo_filename,
    kmers,
    z_thresh,
    pval_thresh
){
  
  for (kmer_idx in seq_along(kmers)) {
    
    kmer = kmers[kmer_idx]
    
    output_path = paste0(
      Result_path,
      "kmer_",
      kmer,
      "/"
    )
    
    ######## Load motif results ################
    
    motif_counts <- read.delim(
      paste0(
        output_path,
        zscore_filename
      )
    )
    
    ######## Filter motifs #####################
    
    selected_motif_idx =
      motif_counts$z_score >= z_thresh &
      motif_counts$p_value < pval_thresh
    
    motif_counts =
      motif_counts[selected_motif_idx, ]
    
    ######## Prepare motifs for sequence logo ####
    
    selected_motifs =
      motif_counts[
        ,
        c("motifs", "total_count")
      ]
    
    n_motifs =
      nrow(selected_motifs)
    
    ######## Plot only if motifs are present ####
    
    if (n_motifs > 0){
      
      # Ignore motif abundance and give equal weight
      selected_motifs =
        selected_motifs$motifs
      
      bits_logo_plot =
        ggseqlogo(
          selected_motifs,
          method = "bits"
        )
      
      probability_logo_plot =
        ggseqlogo(
          selected_motifs,
          method = "prob"
        )
      
      combined_logo_plot =
        gridExtra::grid.arrange(
          bits_logo_plot,
          probability_logo_plot
        )
      
      ######## Save plot #########################
      
      ggsave(
        paste0(
          output_path,
          seq_logo_filename
        ),
        plot = combined_logo_plot,
        width = 7,
        height = 7
      )
    }
  }
}

################################ Paths ################################

motif_analysis_paths <- c(
  "Result/figure/2_Motif_analysis/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/"
)

################## Setting the coverage and mutation thresholds ###################

Coverage_thresholds = c(50)
Diff_thresholds = c(0.03)

####################################################################################

kmers = c(5)
Window_lengths = c(100)

####################################################################################

for (path_idx in seq_along(motif_analysis_paths)) {
  
  for (coverage_idx in seq_along(Coverage_thresholds)) {
    
    for (window_length in Window_lengths) {
      
      motif_path = paste0(
        motif_analysis_paths[path_idx],
        "cov_thresh_",
        Coverage_thresholds[coverage_idx],
        "/sites/win_", window_length, "/"
      )
      
      motif_rnd_bg_path = paste0(
        motif_analysis_paths[path_idx],
        "cov_thresh_",
        Coverage_thresholds[coverage_idx],
        "/sites_rnd/win_", window_length, "/"
      )
      
      create_folder(motif_path)
      
      print(motif_path)
      
      for (threshold_idx in seq_along(Diff_thresholds)) {
        
        print(
          Diff_thresholds[threshold_idx]
        )
        
        ###########################################################
        # files
        ###########################################################
        
        Pos_fasta_file = paste0(
          "Pos_mut_thresh_",
          Diff_thresholds[threshold_idx],
          ".fasta.gz"
        )
        
        Pos_motifs_file = paste0(
          "Pos_motifs_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        ###########################################################
        # 1. Get motifs from FASTA
        ###########################################################
        
        print("Getting motifs from fasta...")
        
        get_motifs_from_fasta(
          fasta_path = paste0(
            motif_path,
            Pos_fasta_file
          ),
          Result_path = motif_path,
          filename = Pos_motifs_file,
          kmers = kmers
        )
        
        
        ###########################################################
        # 2. Enriched motifs using shuffled sequences
        ###########################################################
        
        Pos_randomers_file = paste0(
          "Pos_randomers_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        Pos_zscore_file = paste0(
          "Pos_zscore_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        Pos_seqlogo_file = paste0(
          "Pos_seqlogo_thresh_",
          Diff_thresholds[threshold_idx],
          ".pdf"
        )
        
        N = 10 # 100
        z_thresh = 2
        pval_thresh = 0.05
        
        get_enriched_motifs_by_shuffling(
          motif_path,
          Pos_fasta_file,
          Pos_motifs_file,
          Pos_randomers_file,
          Pos_zscore_file,
          Pos_seqlogo_file,
          kmers,
          N,
          z_thresh,
          pval_thresh
        )
        
        ###########################################################
        # 3. Enriched motifs using random background sequences
        ###########################################################
        
        Pos_rnd_bg_file = paste0(
          "Pos_rnd_bg_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        Pos_zscore_rnd_bg_file = paste0(
          "Pos_zscore_rnd_bg_thresh_",
          Diff_thresholds[threshold_idx],
          ".tsv"
        )
        
        Pos_seqlogo_rnd_bg_file = paste0(
          "Pos_seqlogo_rnd_bg_thresh_",
          Diff_thresholds[threshold_idx],
          ".pdf"
        )
        
        fasta_rnd_bg_region = paste0(
          "MT_",
          Diff_thresholds[threshold_idx],
          "/fasta/"
        )
        
        N = 100 # 10
        z_thresh = 2
        pval_thresh = 0.05
        
        get_enriched_motifs_using_rnd_bg(
          motif_path,
          Pos_motifs_file,
          Pos_rnd_bg_file,
          Pos_zscore_rnd_bg_file,
          Pos_seqlogo_rnd_bg_file,
          motif_rnd_bg_path,
          fasta_rnd_bg_region,
          kmers,
          N,
          z_thresh,
          pval_thresh
        )
      }
    }
  }
}

