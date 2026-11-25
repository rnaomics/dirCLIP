# Published pipeline, scripts 1 to 10

Brushed, parameterised versions of the dirCLIP nanopore analysis scripts (originals by K. Woodward).
Input and output locations and tool paths are variables at the top of each shell script, or command
line arguments of the R scripts 4c, 7, 9 and 10; script 8 lists its input files at the top and reads
`OUT_DIR` from the environment. Statistical settings are those of the original scripts, with the
correction to script 7 (two changes) described at the end.

## Requirements

Guppy (settings FLO-MIN106 and SQK-RNA002 in script 1), minimap2, samtools, bedtools, nanopolish,
NanopolishComp, Nanocompore 1.0.x, R (4.1 or later) with tidyverse and glmGamPoi.

## Reference sequences

Script 2 (mapping and nanopolish eventalign) and script 6 (Nanocompore) must use the same transcript
identifiers and coordinate frame, because the signal deviation sites (channel 1) and the read end
sites (channel 2) are merged in script 10 on transcript and position. The original scripts name
`Mus_musculus.GRCm39.combined.fa` in script 2 and `Mus_musculus.GRCm39.cdna.all.fa` in script 6; the
2026 rerun (`../reanalysis_2026/`) uses one reference per library set for both steps.

## Run order

| step | script | output |
|---|---|---|
| 1 | `1_basecalling.sh` | `merged_all_reads.fastq.gz` (pass and fail reads) |
| 2 | `2_mapping_and_nanopolish.sh` | sorted BAM (minimap2 `-ax map-ont -k14 -uf`; unmapped, reverse strand, secondary and supplementary alignments removed with `-F 2324`) and the nanopolish eventalign table (`--scale-events --print-read-names --samples`) |
| 3 | `3_read_coordinate_extraction.sh` | read start and end coordinates (`bedtools bamtobed`) |
| 4a to 4c | `4a_split_reads_randomly.sh`, `4b_split_nanopolish_ouput.sh`, `4c_split_read_end_coordinate_ouput.R` | optional: three pseudo replicates of the reads, of the eventalign table and of the coordinates |
| 5 | `5_run_nanopolishcomp.sh` | kmer level collapsed eventalign (`NanopolishComp Eventalign_collapse`) |
| 6 | `6_run_nanocompore.sh` | Nanocompore `sampcomp` results (`--min_coverage 5 --comparison_methods GMM --max_invalid_kmers_freq 0.6 --logit`) |
| 7 | `7_nanocompore_collapsed_sites_extraction.R` | significant signal deviation positions, regions merged within 10 nt, volcano plot |
| 8 | `8_read_end_generate_counts_matrix.R` | read end count matrix, `read_end_counts_SRSF3vsUVtsp.tsv` |
| 9 | `9_sig_read_end_extraction_and_collapse.R` | significant read end positions (glmGamPoi adjusted p < 0.05 and `lfc_manual` > 2, the log2 ratio of mean raw counts per condition with one pseudo count each), regions merged within 5 nt, volcano plot |
| 10 | `10_combine_sig_collapsed_read_end_nanocompore.R` | merged site table of both channels, regions of the two channels merged within 5 nt |

Channel 1 (scripts 5 to 7) and channel 2 (scripts 3, 8 and 9) are independent; both feed script 10.

## Example invocations

```bash
# step 1
FAST5_DIR=/data/run/fast5 FASTQ_DIR=/data/run/fastq GUPPY=/opt/guppy/bin/guppy_basecaller bash 1_basecalling.sh

# step 2
REF=/refs/Mus_musculus.GRCm39.cdna.all.fa FASTQ=/data/run/fastq/merged_all_reads.fastq.gz \
  FAST5_DIR=/data/run/fast5 SEQ_SUM=/data/run/sequencing_summary.txt \
  NANO_OUT=/data/run/sample_nano.tsv bash 2_mapping_and_nanopolish.sh

# step 7: SRSF3 (defaults: adjusted p < 0.005, |LOR| >= 1, merge within 10 nt)
Rscript 7_nanocompore_collapsed_sites_extraction.R nanocompore_results.tsv out_dir/
# step 7: HMGA1a (adjusted p < 0.05)
Rscript 7_nanocompore_collapsed_sites_extraction.R nanocompore_results.tsv out_dir/ 0.05 1 10

# step 8: edit the file list at the top, then
OUT_DIR=out_dir Rscript 8_read_end_generate_counts_matrix.R

# step 9
Rscript 9_sig_read_end_extraction_and_collapse.R out_dir/read_end_counts_SRSF3vsUVtsp.tsv out_dir/ UV SRSF3

# step 10
Rscript 10_combine_sig_collapsed_read_end_nanocompore.R \
  out_dir/sig_sites_collapsed_regions_nano.tsv \
  out_dir/sig_sites_collapsed_regions_read_end.tsv \
  out_dir/merged_nano_read_end_sites.tsv
```

## Coordinates and regions

Nanocompore `pos` is the 0 based transcript position of the first base of the 5 nt kmer. Script 7
adds 3, which gives the 1 based position of the kmer centre; this is the `pos` of Supplementary
Table 1. The `GMM_logit_pvalue` column of Nanocompore 1.0.x is already Benjamini Hochberg adjusted
(`../nanocompore_pvalue/`). Scripts 7 and 9 report the middle of each merged span (which may end in
.5) and its window of plus and minus 5; script 10 floors the centre of each merged region.
Supplementary Table 1 reports the floored centre and centre minus 5 to centre plus 5 (top README,
section 6).

## Changes from the original scripts

Reproducibility changes made in June 2026 when the scripts were brushed:

- All absolute paths (Windows and cluster) replaced by variables or arguments.
- Script 3: Snakemake style `{input}` and `{output}` braces, which bash does not expand, changed to `"$input"` and `"$output"`.
- Script 4a: the shuffle step read a file that was never created; the chain of deduplication, shuffle and split now uses one consistent file name.
- Script 4c renamed from `.sh` to `.R` (it contains R code) and given a `write.table` call; the original wrote nothing.
- Scripts 7 and 9: the volcano plots, built but never written, are saved with `ggsave`.

Changes of 10 Oct 2026:

- Script 7: significant positions are merged within `MERGE_GAP` = 10 nt (archived: 5 nt), the value
  that reproduces the signal deviation regions of Supplementary Table 1, and the thresholds are
  optional arguments (defaults 0.005 and 1). Agreed between the coauthors on 5 and 6 Oct 2026.
  Scripts 9 and 10 keep 5 nt, which is what the table reports.
- Script 4c: the subset table is written without a header row (`col.names = FALSE`), because script 8
  reads it with its own column names; with the header the pseudo replicate path failed.
- Script 6: the comment on `--logit` corrected.
