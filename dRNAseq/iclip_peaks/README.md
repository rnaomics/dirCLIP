# iCLIP peaks: liftover, tiers and the depth lookup

Prepares the published SRSF3 iCLIP peaks for comparison with the direct RNA positions, and builds
the position level lookup that carries the dirCLIP read depth. Tools: R2Dtool 1.0.0, UCSC liftOver
with `mm9ToMm39.over.chain.gz`, bedtools 2.31.1. Chromosome names are harmonised to Ensembl style
(`chr1` to `1`, `chrM` to `MT`). Step 4 is a Slurm job: create `logs/` in the submission directory
first, give the partition with `sbatch -p`, and export `CONDA_SH` and `CONDA_ENV` (top README, section 2).

| step | script | what it does |
|---|---|---|
| 1 | `01_prep_r2dtool_input.py --results <nanocompore_results.tsv> --shift <nanocompore_shift_stats.tsv> --out r2dtool_input.tsv` | every row of the SRSF3 Nanocompore results table (2026 rerun), with or without a p value, as the 5 nt kmer footprint (0 based start = pos, end = pos + 5) with payload columns: kmer, p value, LOR, the intensity and dwell metrics, and `test_depth` from `cluster_counts`, summed over the samples whose label contains `_CLIP_` (the test group of the rerun). `test_depth` is 0 where Nanocompore fitted a single cluster; that 0 is a placeholder, not a depth. |
| 2 | `BASE=... INPUT=r2dtool_input.tsv OUT=. bash 02_liftover_and_intersect.sh` (also `R2D`, `GTF`, `ICLIP_BED`, `LIFTOVER`, `CHAIN`) | R2Dtool liftover of the positions to GRCm39 genomic coordinates (`sites_genomic.raw.tsv`, `sites_genomic.sorted.bed`); iCLIP peaks to BED6 (mm9), UCSC liftOver to mm39 (`iclip.mm39.bed`, unmapped peaks in `iclip.mm9_unmapped.bed`), then two variants: whole peak (`iclip.full.sorted.bed`) and 1 nt peak centre (`iclip.centre.sorted.bed`); stranded intersect counting the peaks overlapping each position (`sites_labeled.{full,centre}.bed`). |
| 3 | `OUT=. bash 03_peak_tiers.sh` | confidence tiers by crosslink score: all peaks (370,675), score at or above 15 (38,724 peaks, the top 10.4 percent), score at or above 79 (3,754 peaks, the top 1.0 percent). |
| 4 | `B=<directory of step 2> sbatch 04_build_depth_lookup.slurm` | reformats the lifted positions keeping the transcript position, intersects them with each tier of both variants, and writes `a4a_lookup/a4a_lookup.tsv.gz`: transcript, tx_pos, ref_kmer, our_GMM_p, our_LOR, abs_d_intensity, abs_log2_dwell_fc, test_depth, and the six overlap counts ov_full_all, ov_full_hc15, ov_full_hc79, ov_centre_all, ov_centre_hc15, ov_centre_hc79. |

The lookup of step 4 is the depth source of `../iclip_overlap_published_sites/`, where `test_depth`
of 0 is treated as missing.
