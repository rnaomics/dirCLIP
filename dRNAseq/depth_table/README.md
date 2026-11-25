# Sequencing depth and site yield table

1. `ZIP=/data/fastq.zip OUT=. bash fastq_stats.sh`: streams every `.fastq.gz` member of the zip
   archive (no extraction) and writes, per library, reads, total bases, mean and median length, N50
   and maximum length into `per_sample/` and `depth_table.tsv`. Needs unzip and gawk.
2. `../reanalysis_2026/map_dRNA.slurm` and `gather_counts.sh` give `aligned_read_counts.tsv`, the
   primary alignment counts.
3. `python3 make_dRNA_depth_table.py <TSV_DIR> depth_table.tsv aligned_read_counts.tsv <out_dir> [<T3_DIR>]`
   (`<out_dir>` must exist): joins the two tables with the counts of significant positions and regions
   read from the Supplementary Table 1 TSVs (`<TSV_DIR>`, see `../supplementary_table1_checks/`) and
   writes the supplementary workbook with its notes. With `<T3_DIR>`, the sheets of Supplementary Table 3
   extracted the same way, it adds the HMGA1a site yield. Needs pandas and openpyxl.

This is the 18 Nov 2026 version of the table builder. Against the version of 27 Sep 2026, note 2 states
the basecalls of the FASTQ files used (Guppy 6.1.2 for seven libraries, Dorado 0.2.1 for two) and note 8
the HMGA1a site yield, now a sheet of its own.
