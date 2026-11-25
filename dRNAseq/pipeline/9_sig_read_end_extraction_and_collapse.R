# 9_sig_read_end_extraction_and_collapse.R
# Gamma-Poisson GLM on read-end frequencies between conditions; collapse sites within 5 nt.
# Brushed: paths/conditions are arguments / top-of-file config. Statistics UNCHANGED.
suppressMessages({ library(tidyverse); library(glmGamPoi) })

# ---- inputs: Rscript 9_*.R <counts_matrix.tsv> <output_dir> [ref_cond] [test_cond] ----
args <- commandArgs(trailingOnly = TRUE)
counts_matrix_path <- if (length(args) >= 1) args[1] else "read_end_counts_SRSF3vsUVtsp.tsv"
output_dir         <- if (length(args) >= 2) args[2] else "test_sig_site_generation"
# condition_levels[1] = reference level in the GLM; must match column-name prefixes
condition_levels <- c(if (length(args) >= 3) args[3] else "UV",
                      if (length(args) >= 4) args[4] else "SRSF3")

# 2. load counts, merge Transcript+Start into one ID
counts_raw <- read_tsv(counts_matrix_path, show_col_types = FALSE) %>%
  unite(name, Transcript, Start)
count_cols <- colnames(counts_raw)[colnames(counts_raw) != "name"]

sample_info <- tibble(
  sample = count_cols,
  condition = factor(case_when(
    str_detect(sample, condition_levels[1]) ~ condition_levels[1],
    str_detect(sample, condition_levels[2]) ~ condition_levels[2]
  ), levels = condition_levels)
)
cat("Detected samples and conditions:\n"); print(sample_info)

# 3. filter low-count rows
counts_filtered <- counts_raw %>%
  filter(rowSums(across(all_of(count_cols))) > 10,
         rowSums(across(all_of(count_cols), ~ .x > 0)) >= 2)
cat("Rows after filtering:", nrow(counts_filtered), "\n")

# 4. gamma-Poisson GLM
glm_matrix <- as.matrix(counts_filtered[, count_cols])
rownames(glm_matrix) <- counts_filtered$name
fit <- glm_gp(glm_matrix, design = ~ condition, col_data = sample_info)
contrast_term <- paste0("condition", condition_levels[2])
res <- test_de(fit, contrast = contrast_term) %>% as_tibble() %>%
  mutate(name = rownames(glm_matrix)) %>% relocate(name)

# 5. manual library-normalised LFC (GLM caps extreme LFC at 10; manual better for plotting)
lfc_manual <- counts_filtered %>%
  pivot_longer(cols = all_of(count_cols), names_to = "sample", values_to = "count") %>%
  left_join(sample_info, by = "sample") %>%
  group_by(name, condition) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_count) %>%
  mutate(lfc_manual = log2((.data[[condition_levels[2]]] + 1) /
                           (.data[[condition_levels[1]]] + 1)))

# 6. merge + define significance (adj_pval < 0.05 & lfc_manual > 2)
res_full <- res %>%
  left_join(as_tibble(glm_matrix, rownames = "name"), by = "name") %>%
  left_join(lfc_manual, by = "name") %>%
  separate(name, into = c("tsp", "start"), sep = "_(?=[0-9]+$)") %>%
  mutate(start = as.numeric(start),
         Sig = case_when(adj_pval < 0.05 & lfc_manual > 2 ~ "S", TRUE ~ "NS"))

plot_read_end <- ggplot(res_full %>% filter(lfc_manual > 0),
       aes(x = lfc_manual, y = -log(adj_pval, 10), color = Sig)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(name = "Significance",
                     values = c("S" = "#72C1B4", "NS" = "lightgrey")) +
  geom_hline(yintercept = -log(0.05, 10), colour = "#d7191c", alpha = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 2, colour = "#d7191c", alpha = 0.7, linetype = "dashed") +
  theme_classic() +
  xlab(expression(log[2] * "(fold change)")) + ylab(expression(-log[10](p-value)))

# 7. collapse significant sites within 5 nt
collapsed_regions <- res_full %>%
  filter(Sig != "NS") %>%
  arrange(tsp, start) %>% group_by(tsp) %>%
  mutate(group = cumsum(c(0, diff(start)) > 5)) %>%
  group_by(tsp, group) %>%
  summarise(region_start = min(start), region_end = max(start),
            site_count = n(), max_lfc = max(lfc_manual, na.rm = TRUE),
            min_adj_pval = min(adj_pval, na.rm = TRUE), .groups = "drop") %>%
  mutate(region_middle = (region_start + region_end) / 2,
         window_minus_5 = region_middle - 5, window_plus_5 = region_middle + 5)

# 8. save
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(res_full,                          file.path(output_dir, "GLM_results_sig_marked.tsv"))
write_tsv(collapsed_regions,                 file.path(output_dir, "sig_sites_collapsed_regions_read_end.tsv"))
write_tsv(res_full %>% filter(Sig != "NS"),  file.path(output_dir, "sig_sites_only.tsv"))
ggsave(file.path(output_dir, "read_end_volcano.pdf"), plot_read_end, width = 5, height = 4)
cat("Outputs written to:", output_dir, "\n")
