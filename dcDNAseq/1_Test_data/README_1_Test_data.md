# 1_Test_data

## Purpose

The files in this directory provide small test datasets that allow users to verify that the dcDNAseq workflow runs correctly from input FASTQ files through mutation-site detection and downstream analyses.

These datasets are intended for **pipeline testing and demonstration only** and are not intended to reproduce the biological results presented in the manuscript.

---

## Test FASTQ Files

### SRSF3 dirCLIP dcDNA

```text
SRSF3_dcDNA/
├── rep1_1pct.fastq.gz
├── rep2_1pct.fastq.gz
└── rep3_1pct.fastq.gz
```

Each file contains approximately **1% of reads** randomly sampled from:

```text
SRSF3_DIRCLIP_dcDNA_rep1.fastq.gz
```

Generation commands:

```bash
seqtk sample -s100 SRSF3_DIRCLIP_dcDNA_rep1.fastq.gz 0.01 | gzip > Tools/1_Test_data/SRSF3_dcDNA/rep1_1pct.fastq.gz
seqtk sample -s200 SRSF3_DIRCLIP_dcDNA_rep1.fastq.gz 0.01 | gzip > Tools/1_Test_data/SRSF3_dcDNA/rep2_1pct.fastq.gz
seqtk sample -s300 SRSF3_DIRCLIP_dcDNA_rep1.fastq.gz 0.01 | gzip > Tools/1_Test_data/SRSF3_dcDNA/rep3_1pct.fastq.gz
```

### TotalRNA dcDNA (-UV)

```text
TotalRNA_dcDNA_no_UV/
├── rep1_0.2pct.fastq.gz
├── rep2_0.2pct.fastq.gz
└── rep3_0.2pct.fastq.gz
```

Each file contains approximately **0.2% of reads** randomly sampled from:

```text
SRSF3-EGFP_BAC_Total_RNA_non-UV_dcDNA_rep1.fastq.gz
```

## Running the Test Data

Run the standard workflow beginning from:

```text
3_Genome_Mapping
4_Isoform_Mapping
5_Mutation_profiling
```

## Expected Results

Using the manuscript settings:

```text
Coverage threshold = 50
Mutation-rate difference threshold = 0.03
```

only a small number of positions are expected to satisfy the filtering criteria.

The datasets demonstrate the mechanics of the workflow rather than the biological conclusions of the study.
