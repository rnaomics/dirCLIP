# 8_read_end_generate_counts_matrix.R
# Generate the counts matrix of read-end positions (per Transcript, Start) across replicates
# of two conditions. Brushed: file lists are read from a top-of-file config block (edit these
# to point at your replicate files). Logic UNCHANGED.
suppressMessages(library(tidyverse))

# ---------------- CONFIG: two conditions, each a vector of replicate read-coordinate files ----
# The first list entry is the reference/test condition that defines the position universe.
condition_files <- list(
  SRSF3 = c(
    "read_end_points_SRSF3_rep1.txt",
    "read_end_points_SRSF3_rep2.txt"
  ),
  UV = c(
    "UV_tsp_rep1.txt",
    "UV_tsp_rep2.txt",
    "UV_tsp_rep3.txt"
  )
)
out_dir  <- Sys.getenv("OUT_DIR", unset = ".")
out_file <- "read_end_counts_SRSF3vsUVtsp.tsv"
# ----------------------------------------------------------------------------------------------

# load one replicate -> read-end counts per (Transcript, Start)
load_replicate <- function(filepath, col_name) {
  read_tsv(filepath, col_names = c("Transcript", "Start", "End", "Read_ID", "Strand"),
           show_col_types = FALSE) %>%
    mutate(Start = Start + 1) %>%
    group_by(Transcript) %>% count(Start) %>% ungroup() %>%
    rename(!!col_name := n)
}

# load all replicates for one condition -> wide table
load_condition <- function(filepaths, condition_name) {
  rep_names <- paste0(condition_name, "_rep", seq_along(filepaths))
  map2(filepaths, rep_names, load_replicate) %>%
    reduce(full_join, by = c("Transcript", "Start"))
}

# build full matrix (left_join prioritises positions present in the first condition)
counts_matrix <- imap(condition_files, load_condition) %>%
  reduce(left_join, by = c("Transcript", "Start")) %>%
  replace(is.na(.), 0)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(counts_matrix, file.path(out_dir, out_file))
cat("Wrote", nrow(counts_matrix), "positions to", file.path(out_dir, out_file), "\n")
