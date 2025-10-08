# Set library path
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

# Load libraries
library(ShortRead)
library(stringr)
library(stringi)
library(ggplot2)
library(ggseqlogo)     # Installed via devtools::install_github("omarwagih/ggseqlogo")
library(Biostrings)     # For sequence and motif manipulation
library(tidyverse)
library(writexl)
library(dplyr)

# Set working directory
Project_path <- "/data/groups/anko/mohamed/oxford_nanopore/"
Tools_path <- "Tools/DIRCLIP_analysis/5_Motif_analysis/4_Find_motifs_from_mutation_sites/"
setwd(Project_path)

# Load custom motif functions
# source(paste0(Project_path, Tools_path, "Motif_functions.R"))

## 🧬 Extract Motifs from FASTA

get_motifs_from_fasta <- function(fasta_path, Result_path, filename, kmers) {
  fasta <- readFasta(dirPath = fasta_path)
  fasta_Reads <- data.frame(fasta@sread)[[1]]
  
  for (kmer in kmers) {
    motifs <- find_motifs(kmer, fasta_Reads)
    res_path <- paste0(Result_path, "kmer_", kmer, "/")
    create_folder(res_path)
    
    write.table(motifs, paste0(res_path, filename), sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)
  }
}


##🔎 Identify Most Frequent Motifs

find_motifs <- function(kmer, Reads) {
  if (isEmpty(Reads)) {
    return(data.frame(motifs = "Empty", total_count = 0, prop = 0))
  }
  
  lengths <- nchar(Reads)
  Uniq_lengths <- unique(lengths)
  
  all_kmers <- unlist(lapply(Uniq_lengths, function(L) {
    Reads2 <- Reads[lengths == L]
    n_kmer <- L - kmer + 1
    as.vector(sapply(1:n_kmer, function(x) substr(Reads2, x, kmer + x - 1)))
  }))
  
  df <- data.frame(motifs = all_kmers) %>%
    group_by(motifs) %>%
    summarise(total_count = n(), .groups = 'drop') %>%
    arrange(desc(total_count)) %>%
    mutate(prop = total_count / sum(total_count))
  
  return(df)
}


## 🎲 Shuffle Reads & Compare Motifs
get_motifs_after_shuffling <- function(fasta_path, Result_path, motif_filename,
                                       randomer_filename, kmers, N) {
  fasta <- readFasta(dirPath = fasta_path)
  fasta_Reads <- data.frame(fasta@sread)[[1]]
  
  for (kmer in kmers) {
    motifs <- read.delim(paste0(Result_path, "kmer_", kmer, "/", motif_filename))
    rand_motifs <- get_random_motif_counts(fasta_Reads, motifs, N, kmer)
    
    res_path <- paste0(Result_path, "kmer_", kmer, "/")
    create_folder(res_path)
    
    write.table(rand_motifs, paste0(res_path, randomer_filename), sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


##♻️ Compute Motif Counts After Randomization
get_random_motif_counts <- function(Reads, Top_motifs, N, kmer) {
  Top_motifs <- Top_motifs[,-3]  # Remove proportion column
  Motifs <- Top_motifs
  
  for (kk in 1:N) {
    Reads2 <- stri_rand_shuffle(Reads)
    motifs_rand <- find_motifs(kmer, Reads2)
    
    missing <- setdiff(Motifs$motifs, motifs_rand$motifs)
    if (length(missing) > 0) {
      motifs_rand <- rbind(motifs_rand, data.frame(motifs = missing,
                                                   total_count = 0, prop = 0))
    }
    
    motifs_rand <- motifs_rand[,-3]
    Top_motifs <- merge(Top_motifs, motifs_rand, by = "motifs")
    print(kk)
  }
  
  colnames(Top_motifs)[-1] <- c("T", paste0("C", 1:N))
  Top_motifs <- Top_motifs %>%
    arrange(desc(T)) %>%
    select(-T)
  
  return(Top_motifs)
}


## 📊 Calculate Z-Scores for Motifs
get_motifs_z_score <- function(Result_path, motif_filename,
                               randomer_filename, zscore_filename, kmers) {
  for (kmer in kmers) {
    res_path <- paste0(Result_path, "kmer_", kmer, "/")
    motifs <- read.delim(paste0(res_path, motif_filename))
    randomers <- read.delim(paste0(res_path, randomer_filename))
    
    M_rand <- rowMeans(randomers[,-1])
    SD_rand <- matrixStats::rowSds(as.matrix(randomers[,-1]))
    
    Z_score <- (motifs$total_count - M_rand) / SD_rand
    
    motifs$shuffled_mean <- M_rand
    motifs$shuffled_std <- SD_rand
    motifs$z_score <- Z_score
    
    motifs <- motifs %>% arrange(desc(z_score))
    create_folder(res_path)
    
    write.table(motifs, paste0(res_path, zscore_filename), sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

####🧪 Run Fisher's Exact Test to Assess Motif Enrichment
get_motifs_pvalue <- function(Result_path, zscore_filename, kmers) {
  
  # ⚗️ Function to compute right-tailed Fisher's exact test for each motif
  perform_fisher_test <- function(row) {
    contingency_table <- matrix(row, nrow = 2, byrow = TRUE)
    result <- fisher.test(contingency_table, alternative = "greater")
    return(result$p.value)
  }
  
  for (kmer in kmers) {
    res_path <- paste0(Result_path, "kmer_", kmer, "/")
    motifs <- read.delim(paste0(res_path, zscore_filename))
    
    Sum_total <- sum(motifs$total_count)
    Sum_shuffle <- sum(motifs$shuffled_mean)
    
    df <- motifs[, c("total_count", "shuffled_mean")]
    df$sub_total_count <- Sum_total - df$total_count
    df$sub_shuffle <- Sum_total - df$shuffled_mean
    
    # 🧮 Apply Fisher's exact test row-wise
    motifs$p_value <- apply(df, 1, perform_fisher_test)
    
    motifs <- motifs %>% arrange(p_value)
    
    create_folder(res_path)
    write.table(motifs, paste0(res_path, zscore_filename), sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


### 🎨 Plot Sequence Logos for Statistically Enriched Motifs
make_motifs_sequence_logo_plot <- function(Result_path, zscore_filename, seq_logo_filename,
                                           kmers, z_thresh, pval_thresh) {
  
  for (kmer in kmers) {
    res_path <- paste0(Result_path, "kmer_", kmer, "/")
    motifs <- read.delim(paste0(res_path, zscore_filename))
    
    # 🎯 Filter by z-score and p-value thresholds
    motifs <- motifs %>%
      filter(z_score >= z_thresh & p_value < pval_thresh)
    
    if (nrow(motifs) > 0) {
      Motifs <- motifs$motifs
      
      # 📊 Generate sequence logos using ggseqlogo
      p1 <- ggseqlogo(Motifs, method = 'bits')
      p2 <- ggseqlogo(Motifs, method = 'prob')
      p_combined <- gridExtra::grid.arrange(p1, p2)
      
      ggsave(filename = paste0(res_path, seq_logo_filename),
             plot = p_combined, width = 7, height = 7)
      
      # 🧹 Cleanup session graphics cache (RStudio specific)
      session_path <- "/home/mohamedb/.local/share/rstudio/sessions/active"
      session_name <- list.files(path = session_path)
      unlink(file.path(session_path, session_name, "graphics-r3"), recursive = TRUE)
    }
  }
}

################################################################################################################################################

# 📂 Define base directories and references
paths <- c(
  "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/CLIP_cDNA_merge/"
)

ref_subfolder <- c(
  "ref_P19_WT_TotalRNA_cDNA_no_UV",
  "ref_TotalRNA_cDNA_UV",
  "ref_P19_WT_TotalRNA_cDNA_no_UV",
  "ref_TotalRNA_cDNA_UV"
)

# 📁 Subfolder structures
subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

# 🔢 Threshold parameters
Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
sub_folder_suffix <- ""

# 🧬 Motif extraction parameters
kmers <- c(5)
Lengths <- c(100)
offset <- 0
mut_types <- c("sites", "regions", "regions_centre")

### 🔁 Run Motif Analysis Across Conditions
for (sub1 in subfolders1) {
  for (sub2 in subfolders2) {
    
    Motifs_path <- paste0(paths, sub1, "/", sub2, "/", ref_subfolder, "/")
    
    for (i in seq_along(Motifs_path)) {
      for (coverage in Coverage_thresholds) {
        for (L in Lengths) {
          for (mut_type in mut_types) {
            
            # 📌 Construct the motif sub-directory path
            motif_path <- if (offset == 0) {
              paste0(Motifs_path[i], sub_folder_suffix,
                     "cov_thresh_", coverage, "/win_", L, "/", mut_type, "/")
            } else {
              paste0(Motifs_path[i], sub_folder_suffix,
                     "cov_thresh_", coverage, "/offset_", offset,
                     "_win_", L, "/")
            }
            
            print(motif_path)
            
            for (diff in Diff_thresholds) {
              print(diff)
              
              # 📎 Filenames based on thresholds
              fasta_file     <- paste0("Pos_mut_thresh_", diff, ".fasta.gz")
              motifs_file    <- paste0("Pos_motifs_thresh_", diff, ".tsv")
              randomers_file <- paste0("Pos_randomers_thresh_", diff, ".tsv")
              zscore_file    <- paste0("Pos_zscore_thresh_", diff, ".tsv")
              seqlogo_file   <- paste0("Pos_seqlogo_thresh_", diff, ".pdf")
              
              # ⚙️ Step 1: Extract motifs
              print("Getting motifs from fasta...")
              get_motifs_from_fasta(
                fasta_path = paste0(motif_path, fasta_file),
                Result_path = motif_path,
                filename = motifs_file,
                kmers = kmers
              )
              
              # ⚙️ Step 2: Shuffle reads and get background motifs
              print("Getting random motifs from fasta...")
              set.seed(1)
              get_motifs_after_shuffling(
                fasta_path = paste0(motif_path, fasta_file),
                Result_path = motif_path,
                motif_filename = motifs_file,
                randomer_filename = randomers_file,
                kmers = kmers,
                N = 100
              )
              
              # ⚙️ Step 3: Compute z-score
              print("Calculate the z-score...")
              get_motifs_z_score(
                Result_path = motif_path,
                motif_filename = motifs_file,
                randomer_filename = randomers_file,
                zscore_filename = zscore_file,
                kmers = kmers
              )
              
              # ⚙️ Step 4: Run Fisher's exact test
              print("Calculate the right-tailed Fisher's exact test...")
              get_motifs_pvalue(
                Result_path = motif_path,
                zscore_filename = zscore_file,
                kmers = kmers
              )
              
              # ⚙️ Step 5: Generate logo plots of enriched motifs
              print("Get the sequence logo plot...")
              make_motifs_sequence_logo_plot(
                Result_path = motif_path,
                zscore_filename = zscore_file,
                seq_logo_filename = seqlogo_file,
                kmers = kmers,
                z_thresh = 2,
                pval_thresh = 0.05
              )
            }
          }
        }
      }
    }
  }
}


###################################################################################################################





