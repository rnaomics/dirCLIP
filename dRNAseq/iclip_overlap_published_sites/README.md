# iCLIP overlap of the published SRSF3 site sets

Tests whether the positions called significant in Supplementary Table 1 (the +UV and -UV
Nanocompore sheets) overlap published SRSF3 iCLIP peaks more often than the other tested positions,
crude and adjusted for dirCLIP read depth, and how many high confidence iCLIP peaks the site sets
recover against a null of random tested positions. The window is the 11 nt from 5 nt before to 5 nt
after the kmer centre, in genomic coordinates, on the same strand. Depth comes from the lookup of
`../iclip_peaks/` (step 4), with `test_depth` 0 treated as missing. Confidence intervals of the
odds ratios come from resampling loci (sets of transcripts whose tested spans overlap on the same
strand); the number of resamples and the seed are the constants `B` and `SEED` at the top of
`08_bootstrap.py`. The design was written down before the run (5 Oct 2026) and is held in the
project record.

Python needs numpy, pandas, scipy, statsmodels and openpyxl. `<TSV_DIR>` is the output of
`../supplementary_table1_checks/00_extract_sup_table1.py`. `<out_dir>` is one folder that must exist:
put the three files of step 2 into it first, because steps 7 and 9 read them there, and the local
steps write their tables into it with the prefix `r8`.

| step | where | command | output |
|---|---|---|---|
| 1 | local | `python3 01_prep_r2dtool_input.py <TSV_DIR> kat_positions_r2d_input.tsv` | every row of the two Nanocompore sheets with a finite LOR and p value, as a 1 nt R2Dtool input (start = pos minus 1, end = pos) |
| 2 | cluster | `BASE=... sbatch 02_lift_and_intersect.slurm` (also `R2D`, `GTF`; create `logs/` first) | expects the step 1 file in `$BASE/analysis/R8_published_iclip/` and the peak files of `../iclip_peaks/` in `$BASE/analysis/R3.2_iCLIP_signal/`; writes `r8_positions.tsv.gz` (per position: arm, transcript, pos, LOR, p, genomic position, strand, peak counts per tier), `r8_pairs.tsv.gz` (position to peak pairs with peak scores) and `r8_transcript_versions.tsv` |
| 3 | local | `python3 03_r32_stats.py <out_dir>/r8_positions.tsv.gz a4a_lookup.tsv.gz <out_dir>/r8` | `<out_dir>/r8_r32_table.csv`: per arm, threshold (0.005 and 0.05) and tier, the overlap fractions, the crude odds ratio with its Woolf interval and Fisher exact p value, and the Mantel Haenszel odds ratio over depth strata |
| 4 | cluster | `BASE=... bash 04_offbyone_check.sh` | diagnostic: reading each peak 1 nt wider on its low coordinate side reproduces the percentages of a coauthor's R script, whose BED start lacked the +1; its recorded output is `04_offbyone_check_OUTPUT.txt` |
| 5 | local | `python3 05_recovery_stats.py <out_dir>/r8_positions.tsv.gz <out_dir>/r8_pairs.tsv.gz a4a_lookup.tsv.gz <out_dir>/r8` | `<out_dir>/r8_r34_table.csv`: peaks covered by tested positions, peaks recovered by significant positions, and the null from random tested positions |
| 6 | local | `python3 06_version_sensitivity.py <out_dir>/r8_positions.tsv.gz <out_dir>/r8_transcript_versions.tsv <TSV_DIR>` | diagnostic: how many positions lie on transcripts whose version differs from the GTF, and what dropping them changes |
| 8 | local | `python3 08_bootstrap.py <out_dir>/r8_positions.tsv.gz a4a_lookup.tsv.gz <out_dir>/r8_r32_table.csv <out_dir>/r8_r32_bootstrap.csv` | intervals of the crude and Mantel Haenszel odds ratios by resampling transcripts and by resampling loci; asserts that its point estimates equal those of step 3 |
| 9 | local | `python3 09_independent_checks.py <out_dir> a4a_lookup.tsv.gz` | recomputes the Mantel Haenszel estimates with statsmodels and the peak counts with a pandas join, against the tables of steps 3 and 5 |
| 7 | local | `python3 07_make_supp_table.py <out_dir> <workbook.xlsx> a4a_lookup.tsv.gz <TSV_DIR>` | the supplementary workbook; the numbers taken from the result tables are read, not typed, and the fixed constants (resample counts, strata, tool versions, peak set sizes) are stated in the script |

Steps 3, 5, 8 and 9 print their results; keep the printed output next to the tables.
