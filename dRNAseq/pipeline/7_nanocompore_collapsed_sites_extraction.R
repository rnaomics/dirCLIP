# 7_nanocompore_collapsed_sites_extraction.R
# Extract statistically significant signal deviation sites from the Nanocompore output and
# collapse sites within MERGE_GAP nt (10 nt, the value that reproduces Supplementary Table 1).
# 10 Oct 2026, two corrections to the archived script, agreed between the coauthors:
#   (1) the merge distance is 10 nt, not 5: the archived 5 nt collapse does not reproduce the
#       signal deviation regions of Supplementary Table 1, 10 nt reproduces every one of them
#       (see reanalysis_2026/extract_sig_sites_v2.py and supplementary_table1_checks/);
#   (2) the p value and LOR thresholds are optional arguments (0.005 and 1 for SRSF3, as before;
#       0.05 for HMGA1a, as in the Methods).
# The GMM_logit_pvalue column of Nanocompore 1.0.x is already Benjamini Hochberg adjusted.
# The pos + 3 recentring below applies to Nanocompore output only; the pos column of
# Supplementary Table 1 is already the kmer centre and must not be shifted again.
suppressMessages(library(tidyverse))

# ---- inputs: Rscript 7_*.R <nanocompore_output.tsv> <output_dir> [p_threshold] [lor_threshold] [merge_gap] ----
args <- commandArgs(trailingOnly = TRUE)
nanocompore_output <- if (length(args) >= 1) args[1] else "nanocompore_results.tsv"
output_dir         <- if (length(args) >= 2) args[2] else "test_run/"
P_THRESH   <- if (length(args) >= 3) as.numeric(args[3]) else 0.005   # adjusted GMM_logit_pvalue; 0.05 for HMGA1a
LOR_THRESH <- if (length(args) >= 4) as.numeric(args[4]) else 1       # |Logit_LOR|
MERGE_GAP  <- if (length(args) >= 5) as.numeric(args[5]) else 10      # nt, signal deviation channel

# 2. Import nanocompore output, drop unused columns, recentre position to kmer centre (+3)
nanocompore_file <- read_tsv(nanocompore_output, show_col_types = FALSE) %>%
  dplyr::select(-c(
    genomicPos, strand, GMM_cov_type,
    KS_dwell_pvalue, KS_intensity_pvalue,
    GMM_n_clust, cluster_counts, chr
  )) %>%
  mutate(pos = pos + 3)
nanocompore_file$Logit_LOR <- as.numeric(nanocompore_file$Logit_LOR)  # "NC" -> NA

# 3. extract significant sites
sig_sites <- nanocompore_file %>%
  filter(GMM_logit_pvalue < P_THRESH & abs(Logit_LOR) >= LOR_THRESH)
cat("Significant sites detected:", nrow(sig_sites), "\n")

# 4. annotate significance
nanocompore_file <- nanocompore_file %>%
  mutate(Significance = case_when(
    abs(Logit_LOR) >= LOR_THRESH & GMM_logit_pvalue < P_THRESH ~ "S",
    TRUE ~ "NS"
  ))

# 5. volcano plot
volcano_plot <- ggplot(nanocompore_file %>% na.omit(),
                       aes(x = abs(Logit_LOR), y = -log10(GMM_logit_pvalue),
                           color = Significance)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(name = "Significance",
                     values = c("S" = "#8DA0CB", "NS" = "lightgrey")) +
  geom_hline(yintercept = -log10(P_THRESH), colour = "#d7191c", alpha = 0.7) +
  geom_vline(xintercept = LOR_THRESH, colour = "#d7191c", alpha = 0.7) +
  theme_classic() + xlab("|LOR|") + ylab(expression(-log[10](p-value)))

# 7. collapse significant sites within MERGE_GAP nt
collapsed_regions <- nanocompore_file %>%
  filter(Significance != "NS") %>%
  arrange(ref_id, pos) %>%
  group_by(ref_id) %>%
  mutate(group = cumsum(c(0, diff(pos)) > MERGE_GAP)) %>%
  group_by(ref_id, group) %>%
  summarise(region_start = min(pos), region_end = max(pos),
            site_count = n(), max_lor = max(abs(Logit_LOR), na.rm = TRUE),
            min_adj_pval = min(GMM_logit_pvalue, na.rm = TRUE), .groups = "drop") %>%
  mutate(region_middle = (region_start + region_end) / 2,
         window_minus_5 = region_middle - 5, window_plus_5 = region_middle + 5)

# 8. save outputs
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(nanocompore_file, file.path(output_dir, "nanocompore_full_annotated.tsv"))
write_tsv(sig_sites,        file.path(output_dir, "nanocompore_significant_sites_raw.tsv"))
write_tsv(collapsed_regions,file.path(output_dir, "sig_sites_collapsed_regions_nano.tsv"))
ggsave(file.path(output_dir, "nanocompore_volcano.pdf"), volcano_plot, width = 5, height = 4)
cat("Outputs written to:", output_dir, "\n")
