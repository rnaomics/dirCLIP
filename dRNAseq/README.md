# dirCLIP: direct RNA sequencing analysis code

Code of the direct RNA (nanopore) arm of the dirCLIP manuscript, assembled on 10 Oct 2026 and revised on
18 Nov 2026 (section 5, item 11), for the repository folder `dRNAseq`. It has two layers.

1. `pipeline/`: the published analysis pipeline, scripts 1 to 10, from basecalling to the merged
   site table, as written by the first author and brushed for portability in June 2026. Script 7
   carries the correction agreed between the coauthors on 5 and 6 Oct 2026 (section 5).
2. The 2026 reanalysis and the analyses made for the reviewer response: `reanalysis_2026/` (the
   signal level rerun from the FASTQ and FAST5 of all nine libraries on a Slurm cluster),
   `iclip_peaks/` and `iclip_overlap_published_sites/` (comparison with published SRSF3 iCLIP peaks),
   `supplementary_table1_checks/`, `nanocompore_pvalue/`, `pore_usage/` and `depth_table/`.

The scripts outside `pipeline/` are the versions that produced the reported numbers; they differ from
the copies that ran only by the path brushing and the patches listed in section 5. The `pipeline/`
scripts are brushed copies of the first author's originals (changes listed in `pipeline/README.md`)
with the script 7 correction; they were not rerun in that form. `MANIFEST.tsv` lists every file with
its size and sha256. Each folder has its own README with the commands.

## 1. Folder map

| folder | content | runs on |
|---|---|---|
| `pipeline/` | scripts 1 to 10 of the published pipeline, with the optional pseudo replicate split (4a to 4c) | a GPU machine for Guppy, then any Linux machine with the tools of section 2 |
| `reanalysis_2026/` | Slurm scripts of the 2026 rerun: mapping and read counts, nanopolish eventalign and NanopolishComp collapse, Nanocompore at group level and per replicate, and the site extractor | Slurm cluster |
| `iclip_peaks/` | lifting the published SRSF3 iCLIP peaks from mm9 to GRCm39, the confidence tiers, and the position level lookup that carries read depth | cluster (R2Dtool, bedtools, UCSC liftOver) |
| `iclip_overlap_published_sites/` | overlap of the published SRSF3 site sets (Supplementary Table 1) with the iCLIP peaks: crude and depth adjusted odds ratios, locus bootstrap, peak level recovery, the supplementary workbook | steps 2 and 4 on the cluster, the rest local |
| `supplementary_table1_checks/` | extraction of Supplementary Table 1 to TSV, the region rule check, the cross channel merge test | local |
| `nanocompore_pvalue/` | two checks that `GMM_logit_pvalue` of Nanocompore 1.0.x is Benjamini Hochberg adjusted, with the source evidence and the recorded output | local |
| `pore_usage/` | supplementary figure of pore usage over the runs, from MinKNOW run reports | local |
| `depth_table/` | FASTQ statistics and the supplementary table of sequencing depth and site yield | local |

## 2. Software

Versions are those of the 2026 reanalysis. The scripts in `pipeline/` do not pin versions.

| tool | version | used by |
|---|---|---|
| Guppy | settings FLO-MIN106 and SQK-RNA002 in script 1; the two generations of basecalls are described in section 3 | `pipeline/1` |
| Guppy | 6.1.2+e0556ff, model `2020-09-07_rna_r9.4.1_minion_256_8f8fc47b`, for the FASTQ of seven of the nine libraries (section 3) | basecalls used by `reanalysis_2026/` |
| Dorado | 0.2.1, model `rna002_70bps_hac@v3`, for the FASTQ held of two libraries (section 3) | basecalls used by `reanalysis_2026/` |
| minimap2 | 2.29 (r1283) | `pipeline/2`, `reanalysis_2026/map_dRNA.slurm` |
| samtools | 1.21 | `pipeline/2`, `reanalysis_2026/` |
| bedtools | 2.31.1 | `pipeline/3`, `iclip_peaks/`, `iclip_overlap_published_sites/02` |
| nanopolish | 0.14.0 | `pipeline/2`, `reanalysis_2026/00` |
| NanopolishComp | 0.6.12 | `pipeline/5`, `reanalysis_2026/00` |
| Nanocompore | 1.0.4 (the p value adjustment is the same in releases 1.0.1 to 1.0.4; see `nanocompore_pvalue/`) | `pipeline/6`, `reanalysis_2026/10` and `20` |
| R2Dtool | 1.0.0 | `iclip_peaks/02`, `iclip_overlap_published_sites/02` |
| UCSC liftOver | UCSC utilities binary, with `mm9ToMm39.over.chain.gz` | `iclip_peaks/02` |
| R | 4.1 or later with tidyverse and glmGamPoi for `pipeline/`; base R for `nanocompore_pvalue/` | |
| Python | 3.10 or later with numpy, pandas, scipy, statsmodels, openpyxl and matplotlib | all Python scripts |

Slurm: create `logs/` in the submission directory before `sbatch` (only `submit_all.sh` does it). The
scripts name no partition: give it with `sbatch -p <partition>` or export `SBATCH_PARTITION`. Export
`BASE`, `CONDA_SH` (the `etc/profile.d/conda.sh` of your conda installation) and `CONDA_ENV` (the
environment that holds the tools above; Nanocompore had its own environment for `reanalysis_2026/10`
and `20`); `sbatch` passes exported variables to the job. The `module load samtools/1.21` lines are
those of the cluster used.

## 3. Reference data (not included)

- Ensembl GRCm39 transcriptome FASTA: `Mus_musculus.GRCm39.cdna.all.fa` for the SRSF3 libraries,
  and `Mus_musculus.GRCm39.combined.fa`, the concatenation of the cdna.all and ncrna FASTA files, for
  the HMGA1a libraries and the in vitro polyadenylated transcriptome control. `Mus_musculus.GRCm39.114.gtf`
  (Ensembl release 114) for lifting transcript positions to the genome. Transcripts whose version
  differs between Supplementary Table 1 and this GTF are listed in `r8_transcript_versions.tsv`, written
  by `iclip_overlap_published_sites/02_lift_and_intersect.slurm`; `06_version_sensitivity.py` counts
  them and tests the effect of dropping them.
- SRSF3 iCLIP peaks: `SRSF3_merged.3nt_3nt.peaks.bed.gz`, iCount peak calls on mm9 (BED6, stranded,
  UCSC chromosome names), supplied by the coauthors; UCSC chain file `mm9ToMm39.over.chain.gz`.
- Supplementary Table 1 of the manuscript (`Sup_table1.xlsx`, md5 `69eb33a920618a83470a32e33b6327c0`).
- Raw data: FAST5 and FASTQ of the nine direct RNA libraries named in `reanalysis_2026/map_dRNA.slurm`,
  and the MinKNOW run report of each run (one zip per library in `$BASE/Run_reports/`), which
  `reanalysis_2026/00_eventalign_array.slurm` reads for the sequencing summary and `pore_usage/` reads
  for the duty time tables.
- Basecalls. There are two generations. The published site tables (Supplementary Tables 1 and 3) come
  from the first author's run of `pipeline/`; her basecalling scripts call a Guppy 5 installation
  (release not recorded) with the settings FLO-MIN106 and SQK-RNA002. The FASTQ files deposited with
  the paper were basecalled again from the same FAST5 by M. Mohamed Bahrudeen with Guppy 6.1.2+e0556ff,
  the same settings and the model `2020-09-07_rna_r9.4.1_minion_256_8f8fc47b`. Basecalling was not
  repeated for the 2026 reanalysis: it used his Guppy 6 FASTQ for seven libraries and, for
  SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2 and SRSF3-EGFP_BAC_Total_RNA_UV_dRNA_rep1, his Dorado 0.2.1
  basecalls of October 2023 (`rna002_70bps_hac@v3`). Read identifiers match the MinKNOW records one to
  one, so read counts are unaffected; read lengths and aligned read counts of these two libraries are a
  Dorado measure.

## 4. Run order

- A. Published pipeline: `pipeline/` 1 and 2, then 5, 6 and 7 for the signal deviation channel;
  3 (and 4c for pseudo replicates), then 8 and 9 for the read termination channel; 4a and 4b are
  optional; 10 merges the two channels.
- B. 2026 rerun: `reanalysis_2026/map_dRNA.slurm` then `gather_counts.sh`; `submit_all.sh`
  (eventalign array, then Nanocompore array); `20_nanocompore_perrep.slurm`;
  `extract_sig_sites_v2.py` on each Nanocompore results table.
- C. iCLIP peaks: `iclip_peaks/` 01, 02, 03, 04.
- D. Overlap of the published site sets: `supplementary_table1_checks/00_extract_sup_table1.py`,
  then `iclip_overlap_published_sites/` 01, 02 (cluster), 03, 05, 08, 09, 07; 04 (cluster) and 06
  are diagnostics.
- E. Checks and tables: `supplementary_table1_checks/`, `nanocompore_pvalue/`, `depth_table/`, `pore_usage/`.

## 5. Differences from the archived scripts

1. `pipeline/7_nanocompore_collapsed_sites_extraction.R`: the merge distance is `MERGE_GAP` = 10 nt
   (archived: 5 nt; optional argument 5), and the p value and LOR thresholds are optional arguments 3
   and 4 (defaults 0.005 and 1, the SRSF3 values; HMGA1a used 0.05). The archived 5 nt collapse does not
   reproduce the signal deviation regions of Supplementary Table 1; 10 nt reproduces every one of them
   (header of `reanalysis_2026/extract_sig_sites_v2.py`, and `supplementary_table1_checks/`).
2. `reanalysis_2026/extract_sig_sites_v2.py`: optional fourth argument for the p threshold (default
   0.005, as before); the header paragraph on the earlier extractor and the archived script rewritten
   without internal references, and one word of it corrected ("raw"), because the column is adjusted;
   usage line updated. The region table keeps its column name `min_raw_pval`, which holds the minimum
   adjusted p value of the merged positions.
3. `pipeline/4c_split_read_end_coordinate_ouput.R`: `col.names = FALSE`. The June version wrote a header
   row that script 8 then read as data, so the pseudo replicate path could not run.
4. `pipeline/6_run_nanocompore.sh`: the comment on `--logit` corrected (it selects the logistic
   regression of cluster membership on condition; it does not transform the signal).
5. `reanalysis_2026/gather_counts.sh`: the first line is corrected to `#!/bin/bash`.
6. `iclip_peaks/03_peak_tiers.sh`: new file holding the awk filter that ran inside a larger Slurm
   job on 6 Jul 2026, unchanged.
7. `iclip_overlap_published_sites/`: docstrings and comments no longer carry the internal run labels,
   refer to the design as "the preregistered design of 5 Oct 2026 (project record)", name the
   coauthor's R script neutrally, and in 09 name the result tables as written; 07 and 09 read the
   result tables from the outputs folder given on the command line instead of a `run2/` subfolder;
   07 no longer writes an internal provenance note into the notes sheet; 02 takes `R2D` and `GTF`
   from the environment, with the cluster paths as defaults.
8. `supplementary_table1_checks/00_extract_sup_table1.py`: the workbook path is an argument.
9. Wording in comments across the package: em dashes replaced by punctuation; informal references to
   coauthors replaced by neutral wording; the lifted peak files dated 6 Jul 2026 where a comment said
   "August"; `nanocompore_pvalue/bh_check.py` opens a plain `.tsv` as well as `.tsv.gz`.
10. Markers. Scripts brushed on 10 Oct 2026 carry the header line "Brushed copy for the repository,
    10 Oct 2026". The four Slurm scripts brushed in August 2026 (`reanalysis_2026/00`, `10`, `20` and
    `iclip_peaks/04`) carry a "CLUSTER PATH" comment instead. The `pipeline/` scripts carry the notes
    of the June brushing, listed in `pipeline/README.md`.
11. Revision of 18 Nov 2026, before the public release. The six scripts that activate a conda environment
    (`reanalysis_2026/00`, `10` and `20`, `iclip_peaks/04`, `iclip_overlap_published_sites/02` and `04`)
    take `CONDA_SH` and `CONDA_ENV` from the environment instead of the paths of the cluster used, and
    the Slurm scripts no longer name a partition. The basecalling entries of sections 2 and 3 state the
    two generations of basecalls. `depth_table/make_dRNA_depth_table.py` has new notes 2 and 8 and an
    optional fifth argument for the HMGA1a site yield. Section 6 gained the U1 rule and the position
    convention of Supplementary Tables 3 and 4. One error message of `reanalysis_2026/10` lost its dash.
    No analysis logic changed.

## 6. Conventions

- Nanocompore `pos` is the 0 based transcript position of the first base of the 5 nt kmer. Script 7
  recentres it with pos + 3, which is the 1 based position of the kmer centre. The `pos` column of
  Supplementary Table 1 is already that centre and must not be shifted again.
- `GMM_logit_pvalue`, as written by Nanocompore 1.0.x, is already Benjamini Hochberg adjusted
  across the tested positions (`nanocompore_pvalue/SOURCE_EVIDENCE.md`). Supplementary Table 1 calls
  the column `GMM_adj_pval`.
- Site thresholds. Signal deviation: adjusted p < 0.005 and |LOR| >= 1 for SRSF3, adjusted p < 0.05
  and |LOR| >= 1 for HMGA1a. Read termination: adjusted p < 0.05 and `lfc_manual` > 2, the log2 ratio
  of the mean raw read end counts per condition with one pseudo count each (script 9). Merging: 10 nt
  within the signal deviation channel, 5 nt within the read termination channel, 5 nt between the
  channels (script 10).
- U1 snRNA variants (Supplementary Table 4): the sites are listed on adjusted p < 0.05 alone; the
  |LOR| >= 1 threshold of the other HMGA1a sites was not applied to them.
- Supplementary Tables 3 and 4 give positions as Nanocompore reports them, the 0 based first base of
  the kmer. Supplementary Table 1 gives the kmer centre, which is that position + 3 (first bullet).
- A region is reported as centre minus 5 to centre plus 5. The centre is the floored mean of the first
  and last merged position, and for a combined region the floored mean of the two centres. Scripts 7
  and 9 write the unfloored middle of the merged span and its window; script 10 floors the centre of
  the combined regions. The rule as the table reports it is implemented exactly in
  `reanalysis_2026/extract_sig_sites_v2.py` (signal deviation channel) and
  `supplementary_table1_checks/regions_check.py` (both channels and their merge).
- BED files are 0 based and half open; R2Dtool takes a 0 based start. The iCLIP overlap window is
  the 11 nt from 5 nt before to 5 nt after the kmer centre, on the same strand.

## 7. Supplementary Table 1 against the stated rules

`supplementary_table1_checks/regions_check.py` applies the stated rules (signal deviation positions with adjusted p < 0.005 and |LOR| >= 1 merged within 10 nt; read termination positions with adjusted p < 0.05 and log2 fold change > 2 merged within 5 nt; the two channels merged when their centres are within 5 nt; regions reported as centre minus 5 to centre plus 5) to the position sheets of Supplementary Table 1 and compares the result with the region sheets. Rerun on 10 Oct 2026 with the same output as on 5 Oct 2026.

**+UV sheet.** 742 published rows: 359 signal deviation, 379 read termination and 4 combined regions. The stated rules give 741 regions (358, 378 and 5) and reproduce 740 of the 742 published rows exactly, by transcript, centre and source. The 2 rows not reproduced form 1 pair of neighbouring regions, one from each channel, with centres within 5 nt, which the cross channel rule merges into one combined region:

- ENSMUST00000042235.15: read termination at 1637 and signal deviation at 1638; one combined region at 1637 under the rule

**-UV sheet.** 1,893 published rows: 1,369 signal deviation, 498 read termination and 26 combined regions. Two rows describe the same region (ENSMUST00000042235.15, centre 1637), once as read termination and once as signal deviation, so the sheet holds 1,892 distinct regions. The stated rules give 1,888 regions (1,364, 493 and 31) and reproduce 1,877 of the 1,893 published rows exactly, by transcript, centre and source. Of the 16 rows not reproduced, 14 form 7 pairs of neighbouring regions, one from each channel, with centres within 5 nt, which the cross channel rule merges into one combined region:

- ENSMUST00000042235.15: read termination at 1637 and signal deviation at 1637, the duplicated row; one combined region at 1637 under the rule
- ENSMUST00000058437.14: read termination at 782 and signal deviation at 787; one combined region at 784 under the rule
- ENSMUST00000079419.12: read termination at 615 and signal deviation at 620; one combined region at 617 under the rule
- ENSMUST00000082402.1: read termination at 992 and signal deviation at 996; one combined region at 994 under the rule
- ENSMUST00000082402.1: read termination at 1394 and signal deviation at 1396; one combined region at 1395 under the rule
- ENSMUST00000082402.1: signal deviation at 1420 and read termination at 1424; one combined region at 1422 under the rule
- ENSMUST00000082405.1: signal deviation at 297 and read termination at 302; one combined region at 299 under the rule

The other 2 rows are combined regions whose two components lie more than 5 nt apart, so the rule keeps them separate:

- ENSMUST00000082402.1: combined region at 1046; signal deviation at 1039 and read termination at 1054 under the rule
- ENSMUST00000082402.1: combined region at 1294; signal deviation at 1289 and read termination at 1299 under the rule

The table is kept as submitted, by decision of the coauthors on 6 Oct 2026, and the differences are recorded here. They concern the merging of neighbouring regions only; no position, p value or fold change is affected. Applied strictly, the rules would list 741 regions on the +UV sheet and 1,888 on the -UV sheet.

Assembled by N. Shirokikh, 10 Oct 2026; revised 18 Nov 2026.
