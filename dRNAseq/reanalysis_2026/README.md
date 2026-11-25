# 2026 reanalysis (Slurm cluster)

Rerun of the signal deviation channel for all nine libraries from the FASTQ and the FAST5, made for
the reviewer response. Basecalling was not repeated (top README, sections 2 and 3). Settings are those
of the published pipeline (`../pipeline/`): nanopolish 0.14.0 eventalign with `--scale-events
--print-read-names --samples`, NanopolishComp 0.6.12 `Eventalign_collapse`, Nanocompore 1.0.4
`sampcomp --min_coverage 5 --comparison_methods GMM --max_invalid_kmers_freq 0.6 --logit`. Mapping
used minimap2 2.29 `-ax map-ont -k14 -uf` and samtools 1.21.

Every script takes the working directory from `BASE` and expects `$BASE/ref/` (FASTA, GTF),
`$BASE/fastq/` (one `<library>.fastq.gz` each), `$BASE/fast5/` (one zip archive per library),
`$BASE/Run_reports/` (one MinKNOW run report zip per library, `<library>.zip`, from which step 2 takes
the sequencing summary) and `$BASE/tools/` (minimap2; override with `MINIMAP2`), and writes under
`$BASE/analysis/`. The scripts name no partition (use `sbatch -p` or export `SBATCH_PARTITION`) and
take the conda setup from `CONDA_SH` and `CONDA_ENV` (top README, section 2); the `module load` lines
are those of the cluster used. Create `logs/` in the submission directory before `sbatch`
(`submit_all.sh` does it for the two arrays it submits).

## Libraries and references

| library (FASTQ name) | reference |
|---|---|
| `HMGA1a_DIRCLIP_iv-polyA_dRNA_rep1` | combined |
| `HMGA1a_DIRCLIP_iv-polyA_dRNA_rep2` | combined |
| `HMGA1a_DIRCLIP_iv-polyA_dRNA_rep3` | combined |
| `SRSF3_DIRCLIP_dRNA_rep1` | cdna.all |
| `SRSF3_DIRCLIP_dRNA_rep2` | cdna.all |
| `SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep1` | cdna.all |
| `SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2` | cdna.all |
| `SRSF3-EGFP_BAC_Total_RNA_UV_dRNA_rep1` | cdna.all |
| `WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1` | combined |

`cdna.all` is `Mus_musculus.GRCm39.cdna.all.fa`; `combined` is `Mus_musculus.GRCm39.combined.fa`,
the cdna.all set plus ncrna, used for the HMGA1a libraries and the in vitro polyadenylated control.

## Steps

1. `BASE=... sbatch map_dRNA.slurm` (array 0 to 8): minimap2 to the library's reference, sorted and
   indexed BAM, and per library counts of all records and of primary alignments
   (`samtools view -c -F 0x904`) in `analysis/mapping/bam/<library>.count.tsv`. Then
   `BASE=... bash gather_counts.sh` writes `analysis/mapping/aligned_read_counts.tsv`, an input of
   `../depth_table/`.
2. `BASE=... bash submit_all.sh` (it changes to its own directory, so it must sit next to the two
   array scripts): submits `00_eventalign_array.slurm` (array 0 to 8: FAST5 extraction, sequencing
   summary from the run report, nanopolish index, eventalign, NanopolishComp collapse; a library with an
   existing `out_eventalign_collapse.tsv` is skipped) and, after it, `10_nanocompore_array.slurm`
   (array 0 to 1). Task 0: SRSF3, control = SRSF3 transcriptome without UV (two replicates), test =
   SRSF3 dirCLIP (two replicates), reference cdna.all. Task 1: HMGA1a, control = wild type in vitro
   polyadenylated transcriptome (one replicate), test = HMGA1a dirCLIP (three replicates), reference
   combined. Output: `analysis/nanocompore/<comparison>/<comparison>nanocompore_results.tsv` and
   `..._shift_stats.tsv`.
3. `BASE=... sbatch 20_nanocompore_perrep.slurm` (array 0 to 4): the same comparisons with one test
   replicate at a time against the same control group, for the replicate reproducibility analysis.
4. `python3 extract_sig_sites_v2.py <nanocompore_results.tsv> <out_dir> [prefix] [p_threshold]`
   (create `<out_dir>` first; `prefix` defaults to the input file stem): recentres positions (pos + 3),
   keeps adjusted p below the threshold (default 0.005; 0.05 for HMGA1a) and |LOR| >= 1, merges
   positions within 10 nt on a transcript, reports each region as centre minus 5 to centre plus 5 with
   the floored centre, and also keeps the merged span. Outputs `<prefix>_v2_significant_sites_raw.tsv`
   and `<prefix>_v2_sig_sites_collapsed_regions.tsv`; in the region table the column `min_raw_pval`
   holds the minimum adjusted p value of the merged positions (the name predates the finding that the
   column is adjusted). The header of the script records why the merge distance is 10 nt and not the
   5 nt of the archived R script.
