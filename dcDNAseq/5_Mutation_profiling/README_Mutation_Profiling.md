# 🧬 Mutation Profiling Pipeline

This pipeline performs per-base mutation profiling from Oxford Nanopore sequencing alignments using the `perbase` tool.

## Requirements

- Python 3.7+
- perbase
- pandas
- gzip

---

## Working Directory

```bash
cd Tool_submission/5_Mutation_profiling
```

---

# 1. Run Per-base Profiling

Generate per-base nucleotide counts from BAM files.

```bash
python mutation_profiling_pipeline.py perbase \
--bam_dirs \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--output_dirs \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--genome ../genome/mm9/mm9.fa \
--max_depth 500000 \
--min_base_quality 8 \
--threads 2 \
--cov_thresh 10 \
--mut_thresh 0 \
--min_read_quality 8
```

---

# 2. Calculate Mutation Rates

## Including Deletions

```bash
python mutation_profiling_pipeline.py mutation \
--output_dirs \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--cov_thresh 10 \
--mut_thresh 0 \
--min_base_quality 8 \
--min_read_quality 8 \
--incdel
```

This creates mutation rates where deletions are included in both coverage and mutation counts.

## Excluding Deletions

```bash
python mutation_profiling_pipeline.py mutation \
--output_dirs \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--cov_thresh 10 \
--mut_thresh 0 \
--min_base_quality 8 \
--min_read_quality 8
```

This creates mutation rates based only on nucleotide substitutions.

---

# 3. Run Complete Workflow (Genome-Based Analysis)

## Complete workflow including deletions

```bash
python mutation_profiling_pipeline.py all \
--bam_dirs \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--output_dirs \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--genome ../genome/mm9/mm9.fa \
--max_depth 500000 \
--min_base_quality 8 \
--threads 2 \
--cov_thresh 10 \
--mut_thresh 0 \
--min_read_quality 8 \
--incdel
```

## Complete workflow excluding deletions

```bash
python mutation_profiling_pipeline.py all \
--bam_dirs \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--output_dirs \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--genome ../genome/mm9/mm9.fa \
--max_depth 500000 \
--min_base_quality 8 \
--threads 2 \
--cov_thresh 10 \
--mut_thresh 0 \
--min_read_quality 8
```

---
# 4. Run Complete Workflow (Isoform-Specific Analysis)

## Complete workflow including deletions

```bash
python mutation_profiling_pipeline.py all \
--bam_dirs \
../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--output_dirs \
../Result/5_Mutation_profiling/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--genome ../genome/mm9/Annotation/gencode.vM1.annotation_transcripts_with_chr.fa \
--max_depth 500000 \
--min_base_quality 8 \
--threads 2 \
--cov_thresh 10 \
--mut_thresh 0 \
--min_read_quality 8 \
--incdel
```

## Complete workflow excluding deletions

```bash
python mutation_profiling_pipeline.py all \
--bam_dirs \
../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--output_dirs \
../Result/5_Mutation_profiling/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
../Result/5_Mutation_profiling/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--genome ../genome/mm9/Annotation/gencode.vM1.annotation_transcripts_with_chr.fa \
--max_depth 500000 \
--min_base_quality 8 \
--threads 2 \
--cov_thresh 10 \
--mut_thresh 0 \
--min_read_quality 8
```

---
## Output Structure

For:

```bash
--min_read_quality 8
--min_base_quality 8
```

Results are stored in:

```text
rep1_1pct/
└── q_8_Q_8/
    ├── mutation.txt.gz
    ├── Coverage_10.tsv.gz
    ├── 1_No_del/
    │   └── Mutation_rate_cov_10_mut_0.tsv.gz
    └── 2_Inc_del/
        └── Mutation_rate_cov_10_mut_0.tsv.gz
```

---

## Parameter Description

| Parameter | Description |
|-----------|-------------|
| --bam_dirs | Directory containing BAM files |
| --output_dirs | Output directory |
| --genome | Reference genome FASTA |
| --max_depth | Maximum depth used by perbase |
| --min_base_quality | Minimum base quality threshold |
| --min_read_quality | Determines BAM file (filt_qX.bam) |
| --cov_thresh | Minimum coverage threshold |
| --mut_thresh | Minimum mutation count threshold |
| --threads | Number of threads |
| --incdel | Include deletions in mutation-rate calculation |

---

