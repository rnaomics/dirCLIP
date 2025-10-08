### Set working directory and library path
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### Load required libraries
library(Rsubread)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(S4Vectors)
library(matrixStats)

################################################################################
## Function: write_wiggle2
## Purpose: Write a dataframe in Wiggle format, used for genome browser tracks
################################################################################
write_wiggle2 <- function(df, path, file) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  file_path <- paste0(path, file)
  con <- file(file_path, "w")
  writeLines("track type=wiggle_0", con)
  
  for (chrom in unique(df$chromosome)) {
    writeLines(paste("variableStep chrom=", chrom, sep=""), con)
    df_sub <- df[df$chromosome == chrom, ]
    line_data <- paste(df_sub$start, df_sub$value, sep="\t")
    writeLines(line_data, con)
  }
  
  close(con)
  cat("✔ Wiggle file successfully written to:", file_path, "\n")
}

################################################################################
## Function: add_mut_ratio
## Purpose: Integrate mutation ratios from multiple replicates and compute mean ± SEM
################################################################################
add_mut_ratio <- function(Mut_path, Mut_sample, Mut_bed_filename, cov_thresh, df) {
  for (i in seq_along(Mut_path)) {
    mut <- read.delim(paste0(Mut_path[i], Mut_bed_filename), header = FALSE)
    colnames(mut) <- c("chr", "pos", "cov", "Mut_ratio")
    mut <- mut[mut$cov >= cov_thresh, ]
    rownames(mut) <- paste0(mut$chr, ":", mut$pos)
    mut <- mut[, c("chr", "pos", "Mut_ratio")]
    
    sample_name <- paste0(Mut_sample, "_rep", i)
    colnames(mut)[3] <- sample_name
    
    if (is.null(df)) {
      df <- mut
    } else {
      shared <- intersect(rownames(df), rownames(mut))
      df <- df[shared, ]
      mut <- mut[shared, ]
      df[, sample_name] <- mut[, 3]
    }
  }
  
  reps <- colnames(df) %in% paste0(Mut_sample, "_rep", seq_along(Mut_path))
  df[, paste0(Mut_sample, "_mean")] <- rowMeans(df[, reps])
  df[, paste0(Mut_sample, "_sem")] <- rowSds(as.matrix(df[, reps])) / sqrt(sum(reps))
  return(df)
}

################################################################################
## Function: calculate_p_value
## Purpose: Compute p-values using either t-test or Wilcoxon test
################################################################################
calculate_p_value <- function(x, y, test_type) {
  noise <- runif(length(x), -1e-6, 1e-6)
  x <- x + noise
  y <- y + noise
  
  if (test_type == "t-test") {
    return(t.test(x, y, var.equal = FALSE, paired = FALSE, alternative = "greater")$p.value)
  } else if (test_type == "wilcox-test") {
    return(wilcox.test(x, y, paired = FALSE, alternative = "greater")$p.value)
  }
}

################################################################################
## Function: remove_rand_chr
## Purpose: Filter out entries from random chromosomes
################################################################################
remove_rand_chr <- function(Sites) {
  rand_ids <- grep("random", Sites$V1)
  if (length(rand_ids) > 0) Sites <- Sites[-rand_ids, ]
  return(Sites)
}

################################################################################
## Function: get_subsamples_paths
## Purpose: Construct paths to all subsample folders
################################################################################
get_subsamples_paths <- function(paths, sub_folder) {
  lapply(paths, function(p) paste0(p, "subsample", 1:3, "/", sub_folder, "/"))
}

################################################################################
## Function: get_mut_diff
## Purpose: Calculate differential mutation rates between target and control
################################################################################
get_mut_diff <- function(Target_path, Target_sample, Mut_bed_filename, cov_thresh, Res_path, Res_sub_path, sub_folder_suffix) {
  df <- add_mut_ratio(Target_path, Target_sample, Mut_bed_filename, cov_thresh, NULL)
  df <- add_mut_ratio(Control_path, Control_sample, Mut_bed_filename, cov_thresh, df)
  
  # Compute difference metrics
  df$diff_mean <- df[, paste0(Target_sample, "_mean")] - df[, paste0(Control_sample, "_mean")]
  df$diff_sem <- sqrt(df[, paste0(Target_sample, "_sem")]^2 + df[, paste0(Control_sample, "_sem")]^2)
  
  # p-value computations
  ratio_matrix <- as.matrix(df[, grep("rep", names(df))])
  target_cols <- names(df) %in% paste0(Target_sample, "_rep", seq_along(Target_path))
  control_cols <- names(df) %in% paste0(Control_sample, "_rep", seq_along(Control_path))
  
  df$p_value <- apply(ratio_matrix, 1, function(r) calculate_p_value(r[target_cols], r[control_cols], "t-test"))
  df$wilcox_p_value <- apply(ratio_matrix, 1, function(r) calculate_p_value(r[target_cols], r[control_cols], "wilcox-test"))
  
  # Save result
  out_path <- paste0(Res_path, Res_sub_path, "/", sub_folder_suffix, "cov_thresh_", cov_thresh, "/")
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  write.table(df, paste0(out_path, "significant_mutsites.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  return(df)
}

################################################################################
## Function: filt_mutsites_by_diff_thresh
## Purpose: Filter mutation sites by threshold and export in Wiggle format
################################################################################
filt_mutsites_by_diff_thresh <- function(cov_thresh, Diff_thresholds, wiggle_out, Res_path, Res_sub_path, sub_folder_suffix, filt_by) {
  path <- paste0(Res_path, Res_sub_path, "/", sub_folder_suffix, "cov_thresh_", cov_thresh, "/")
  df <- read.delim(paste0(path, "significant_mutsites.tsv"))
  
  # Add FDR column
  df$FDR <- p.adjust(df$p_value, method = "BH")
  write.table(df, paste0(path, "significant_mutsites.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  pval_pass <- df$p_value < 0.05
  fdr_pass <- df$FDR < 0.05
  
  for (thresh in Diff_thresholds) {
    diff_pass <- df$diff_mean >= thresh
    combined_filter <- if (filt_by == "pvalue") (pval_pass & diff_pass) else (fdr_pass & diff_pass)
    
    filtered_df <- df[combined_filter, ]
    pos_diff <- filtered_df[, c("chr", "pos", "diff_mean")]
    names(pos_diff) <- c("chromosome", "start", "value")
    
    output_path <- paste0(path, filt_by, "/")
    if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
    
    write.table(pos_diff, paste0(output_path, "Pos_mut_thresh_", thresh, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    if (wiggle_out) {
      write_wiggle2(pos_diff, output_path, paste0("Pos_mut_thresh_", thresh, ".wig"))
    }
  }
}

########################################################################################################
## Parameters: Define sample comparisons and thresholds
Target_paths_All <- c(
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/"
)

Target_samples <- c("Hmga1a_CLIP_cDNA", "Hmga1b_CLIP_cDNA")

Control_paths_All <- rep(
  "Results/Nanopore/5_UV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/TotalRNA_cDNA_UV_merge/",
  2
)

Control_samples <- rep("TotalRNA_cDNA_UV", 2)

Res_paths_All <- c(
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/"
)

Res_sub_paths <- paste0("ref_", Control_samples)

subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

## Analysis settings
wiggle_out <- TRUE                # Export Wiggle files for genome browsers
Coverage_thresholds <- c(50, 100, 200, 300, 400, 500)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
filt_by <- "FDR"                 # Filtering method: "pvalue" or "FDR"

## Mutation data filenames
Mut_bed_filename <- "Mutation_rate_cov_10_mut_0.tsv.gz"
All_mut_wig_filename <- "All_mutations.wig"
sub_folder_suffix <- ""

## Main execution loop
for (subfolder1 in subfolders1) {
  for (subfolder2 in subfolders2) {
    Target_paths <- paste0(Target_paths_All, subfolder1, "/")
    Control_paths <- paste0(Control_paths_All, subfolder1, "/")
    Res_paths <- paste0(Res_paths_All, subfolder1, "/", subfolder2, "/")
    
    ## Expand to subsample paths
    Target_paths <- get_subsamples_paths(Target_paths, subfolder2)
    Control_paths <- get_subsamples_paths(Control_paths, subfolder2)
    
    for (i in seq_along(Target_paths)) {
      Target_path <- Target_paths[[i]]
      Target_sample <- Target_samples[i]
      
      Control_path <- Control_paths[[i]]
      Control_sample <- Control_samples[i]
      
      Res_path <- Res_paths[i]
      Res_sub_path <- Res_sub_paths[i]
      
      message("Analyzing comparison: ", Target_sample, " vs ", Control_sample)
      
      for (cov_thresh in Coverage_thresholds) {
        message("Coverage Threshold: ", cov_thresh)
        
        df_mut_ratio <- get_mut_diff(
          Target_path, Target_sample,
          Mut_bed_filename, cov_thresh,
          Res_path, Res_sub_path, sub_folder_suffix
        )
        
        filt_mutsites_by_diff_thresh(
          cov_thresh, Diff_thresholds, wiggle_out,
          Res_path, Res_sub_path, sub_folder_suffix, filt_by
        )
      }
    }
  }
}


