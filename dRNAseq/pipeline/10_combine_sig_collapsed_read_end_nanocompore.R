# 10_combine_sig_collapsed_read_end_nanocompore.R
# Merge significant collapsed sites from the two detection channels (signal-deviation +
# read-end termination) into one site set, collapsing across channels within 5 nt.
# Brushed: paths as arguments / config. Logic UNCHANGED.
suppressMessages(library(tidyverse))

# ---- inputs: Rscript 10_*.R <nano_collapsed.tsv> <readend_collapsed.tsv> <out.tsv> ----
args <- commandArgs(trailingOnly = TRUE)
nano_path    <- if (length(args) >= 1) args[1] else "sig_sites_collapsed_regions_nano.tsv"
readend_path <- if (length(args) >= 2) args[2] else "sig_sites_collapsed_regions_read_end.tsv"
output_path  <- if (length(args) >= 3) args[3] else "merged_nano_read_end_sites.tsv"

shared_cols <- c("ref_id", "region_middle", "minus_5", "plus_5",
                 "site_count", "source", "max_LOR_or_LFC", "min_p_value")

# load both channels, harmonise column names, tag source
combined <- bind_rows(
  read_tsv(nano_path, show_col_types = FALSE) %>%
    rename(minus_5 = window_minus_5, plus_5 = window_plus_5,
           max_LOR_or_LFC = max_lor, min_p_value = min_adj_pval) %>%
    mutate(source = "signal_dev") %>% select(all_of(shared_cols)),
  read_tsv(readend_path, show_col_types = FALSE) %>%
    rename(ref_id = tsp, minus_5 = window_minus_5, plus_5 = window_plus_5,
           max_LOR_or_LFC = max_lfc, min_p_value = min_adj_pval) %>%
    mutate(source = "read_ends") %>% select(all_of(shared_cols))
)

# collapse across channels within 5 nt
collapsed_merged <- combined %>%
  arrange(ref_id, region_middle) %>%
  group_by(ref_id) %>%
  mutate(group = cumsum(c(0, diff(region_middle)) > 5)) %>%
  group_by(ref_id, group) %>%
  summarise(merged_start = min(minus_5), merged_end = max(plus_5),
            region_center = floor((merged_start + merged_end) / 2),
            sources = paste(unique(source), collapse = "+"),
            total_sites = sum(site_count),
            max_LOR_or_LFC = max(max_LOR_or_LFC, na.rm = TRUE),
            min_p_value = min(min_p_value, na.rm = TRUE), .groups = "drop")

write_tsv(collapsed_merged, output_path)
cat("Written:", nrow(collapsed_merged), "collapsed regions to", output_path, "\n")
