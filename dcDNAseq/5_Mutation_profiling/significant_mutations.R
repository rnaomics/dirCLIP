#!/usr/bin/env Rscript
#### Set library path #########
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

# Significant Mutation Site Identification Pipeline
# Refactored for reproducibility and command-line execution

suppressPackageStartupMessages({
  library(optparse)
  library(matrixStats)
  library(dplyr)
})

write_wiggle <- function(df, outfile) {
  con <- file(outfile, 'w')
  writeLines('track type=wiggle_0', con)

  chroms <- unique(df$chromosome)

  for (chr in chroms) {
    writeLines(paste0('variableStep chrom=', chr), con)
    chr_df <- df[df$chromosome == chr, ]
    lines <- paste(chr_df$start, chr_df$value, sep='\t')
    writeLines(lines, con)
  }

  close(con)
}

calculate_p_value <- function(x, y, test_type='t-test') {
  set.seed(1) 
  if (test_type == 't-test') {
    noise <- runif(length(x), -1e-6, 1e-6)
    x <- x + noise
    y <- y + noise

    return(t.test(
      x,
      y,
      alternative='greater'
    )$p.value)
  }

  if (test_type == 'wilcox-test') {
    return(wilcox.test(
      x,
      y,
      alternative='greater'
    )$p.value)
  }
}

add_mut_ratio <- function(paths,
                          sample_name,
                          mutation_file,
                          cov_thresh,
                          df=NULL) {

  for (i in seq_along(paths)) {

    mut <- read.delim(file.path(paths[i], mutation_file), header=FALSE)

    colnames(mut) <- c('chr', 'pos', 'cov', 'Mut_ratio')

    mut <- mut[mut$cov >= cov_thresh, ]

    rownames(mut) <- paste0(mut$chr, ':', mut$pos)

    mut <- mut[, c('chr','pos','Mut_ratio')]

    rep_name <- paste0(sample_name, '_rep', i)
    colnames(mut)[3] <- rep_name

    if (is.null(df)) {
      df <- mut
    } else {

      keep1 <- rownames(mut) %in% rownames(df)
      mut <- mut[keep1, ]

      keep2 <- rownames(df) %in% rownames(mut)
      df <- df[keep2, ]

      df[[rep_name]] <- mut[,3]
    }
  }

  ids <- colnames(df) %in% paste0(sample_name, '_rep', seq_along(paths))

  df[[paste0(sample_name, '_mean')]] <- rowMeans(df[, ids])

  df[[paste0(sample_name, '_sem')]] <- rowSds(as.matrix(df[, ids])) / sqrt(sum(ids))

  return(df)
}

get_mut_diff <- function(target_paths,
                         target_name,
                         control_paths,
                         control_name,
                         mutation_file,
                         cov_thresh,
                         result_dir) {

  df <- NULL
  print(cov_thresh)
  df <- add_mut_ratio(target_paths, target_name, mutation_file, cov_thresh, df)

  df <- add_mut_ratio(control_paths, control_name, mutation_file, cov_thresh, df)

  df$diff_mean <-
    df[[paste0(target_name, '_mean')]] -
    df[[paste0(control_name, '_mean')]]

  df$diff_sem <- sqrt(
    df[[paste0(target_name, '_sem')]]^2 +
    df[[paste0(control_name, '_sem')]]^2
  )

  mut_matrix <- as.matrix(df[, grep('rep', colnames(df))])

  target_cols <- colnames(mut_matrix) %in%
    paste0(target_name, '_rep', seq_along(target_paths))

  control_cols <- colnames(mut_matrix) %in%
    paste0(control_name, '_rep', seq_along(control_paths))

  df$p_value <- apply(mut_matrix, 1, function(row) {
    calculate_p_value(row[target_cols], row[control_cols], 't-test')
  })

  df$FDR <- p.adjust(df$p_value, method='BH')
  
  out_dir = file.path(result_dir, paste0("cov_thresh_", cov_thresh))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  write.table(
    df,
    file.path(out_dir, 'significant_mutsites.tsv'),
    sep='\t',
    quote=FALSE,
    row.names=FALSE
  )

  df
}

filter_sites <- function(
    result_dir,
    cov_thresh,
    diff_threshold,
    wiggle_out = FALSE) {
  
  out_dir <- file.path(result_dir, paste0("cov_thresh_", cov_thresh))
  
  df <- read.delim(file.path(out_dir, "significant_mutsites.tsv"))
  
  keep <- (df$FDR < 0.05 & df$diff_mean >= diff_threshold)
  
  df <- df[keep, c("chr", "pos", "diff_mean")]
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  output_file <- file.path(out_dir, paste0("Pos_mut_thresh_", diff_threshold, ".tsv"))
  
  write.table(
    df,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  if (wiggle_out) {
    
    wig_df <- data.frame(chromosome = df$chr, start = df$pos, value = df$diff_mean)
    
    write_wiggle(
      wig_df,
      file.path(out_dir, paste0("Pos_mut_thresh_", diff_threshold, ".wig"))
    )
  }
}

group_sites_by_region <- function(
    result_dir,
    cov_thresh,
    diff_threshold,
    group_thresh = 5){
  
  input_dir <- file.path(result_dir, paste0("cov_thresh_", cov_thresh))
  
  sites_file <- file.path(input_dir, paste0("Pos_mut_thresh_", diff_threshold, ".tsv"))
  
  region_file <- file.path(input_dir, paste0("Pos_mut_thresh_", diff_threshold, "_region.tsv"))
  print(sites_file)
  sites <- read.delim(sites_file, header = FALSE)
  colnames(sites) = c("chr", "pos", "diff_mean")
  
  sites <- sites[order(sites$chr, sites$pos),]
  
  Sites_grouped <- sites %>%
    dplyr::select(chr, pos) %>% 
    mutate(group = cumsum(diff(c(-Inf, pos)) >= group_thresh | chr != lag(chr, default = first(chr)))) %>%
    group_by(chr, group) %>%
    summarise(start = min(pos),
              end = max(pos),
              range = end - start + 1, .groups = "drop") %>%
    select(chr, start, end, range)
  
  write.table(
    Sites_grouped,
    file = region_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  cat("Saved grouped regions:", region_file, "\n")
}

option_list <- list(
  make_option('--target_paths', type='character'),
  make_option('--control_paths', type='character'),
  make_option('--target_name', type='character'),
  make_option('--control_name', type='character'),
  make_option('--mutation_file', type='character', default='Mutation_rate_cov_10_mut_0.tsv.gz'),
  # make_option('--coverage_threshold', type='integer', default=50),
  # make_option('--diff_threshold', type='double', default=0.03),
  make_option('--coverage_threshold', type='character', default="50"),
  make_option('--diff_threshold', type='character', default="0.03"),
  make_option('--group_threshold', type='integer', default=5),
  make_option('--min_read_quality', type='integer', default=8),
  make_option('--min_base_quality', type='integer', default=8),
  make_option('--result_dir', type='character'),
  make_option('--wiggle', action='store_true', default=FALSE),
  make_option('--incdel', action='store_true', default=FALSE)
)

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

### Create subfolders ############
subfolder1 <- paste0("q_", opt$min_read_quality, "_Q_",  opt$min_base_quality,  "/")

if (opt$incdel) {
  
  subfolder2 <- "2_Inc_del/"
  
} else {
  
  subfolder2 <- "1_No_del/"
}

#################################
if (!is.null(opt$target_paths)) {
  
  target_paths <- strsplit(opt$target_paths, ',')[[1]]
  control_paths <- strsplit(opt$control_paths, ',')[[1]]
  
  target_paths <- paste0(target_paths, subfolder1, subfolder2)
  control_paths <- paste0(control_paths, subfolder1, subfolder2)
  result_dir = paste0(opt$result_dir, subfolder1, subfolder2)
  
  coverage_threshold <- as.integer(strsplit(opt$coverage_threshold, ',')[[1]])
  diff_threshold <- as.double(strsplit(opt$diff_threshold, ',')[[1]])
  
  for (cov_thresh in coverage_threshold){
    df <- get_mut_diff(
      target_paths,
      opt$target_name,
      control_paths,
      opt$control_name,
      opt$mutation_file,
      cov_thresh,
      result_dir
    )
    
    for (diff_thresh in diff_threshold){
      filter_sites(result_dir, cov_thresh, diff_thresh, opt$wiggle)
      
      group_sites_by_region(result_dir, cov_thresh, diff_thresh, opt$group_threshold)
      
    }
  }
}
