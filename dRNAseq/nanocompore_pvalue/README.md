# Is `GMM_logit_pvalue` raw or adjusted?

In Nanocompore 1.0.x the `GMM_logit_pvalue` column of the results table is the Benjamini Hochberg
adjusted p value of the logit test, adjusted in place across all tested positions. `SOURCE_EVIDENCE.md`
gives the lines of the 1.0.4 source (the logit test in TxComp.py with one pseudo read added to each
of the four condition by cluster cells, the adjustment in SampCompDB.py, and the report writer), and
why the development branch reads differently; the adjustment is the same in releases 1.0.1 to 1.0.4.
Two checks:

- `Rscript check_nanocompore_p_is_BH.R <nanocompore_results.tsv or .tsv.gz>` (base R, any 1.0.x
  results table): recomputes the raw logit p value from `cluster_counts` with the pseudo count,
  adjusts it, and compares both with the column. Its output on the 2026 SRSF3 rerun is
  `check_nanocompore_p_is_BH_OUTPUT.txt`.
- `python3 bh_check.py <nanocompore_results.tsv or .tsv.gz>` (numpy): the same question written
  from scratch, and a test of whether the Benjamini Hochberg fixed point argument alone could tell
  raw from adjusted (it cannot). It is written for the SRSF3 comparison: it takes a sample as control
  when its label contains `noUV` and as test when it contains `CLIP`, and skips any other row.

Consequence for the pipeline: the threshold of 0.005 in `../pipeline/7` and in
`../reanalysis_2026/extract_sig_sites_v2.py` is applied to an adjusted p value, and Supplementary
Table 1 is right to call the column `GMM_adj_pval`.
