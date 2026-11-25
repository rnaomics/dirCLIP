# Read Length Statistics

## Purpose

This utility script calculates read length statistics from aligned BAM files.

The script generates:

- Mean read length
- Median read length
- Number of mapped reads

## Requirements

### Software

- samtools
- Python 3

### Python Packages

```bash
pip install pandas openpyxl
```

## Input

Sorted BAM files:

```text
Result/3_Genome_Mapping/*/sorted.bam
```

## Run

```bash
python read_length_statistics.py
```

## Workflow

```text
sorted.bam
    ↓
samtools view
    ↓
Read length distribution
    ↓
Mean read length
Median read length
Read count
    ↓
Read_length.xlsx
```

## Output Files

For each BAM directory:

```text
read_length.txt
Read_length.xlsx
```

### Read_length.xlsx

Contains:

```text
Mean_RL
Median_RL
N_reads
```

## Notes

Reads are extracted using:

```bash
samtools view -F 2308
```

Secondary, supplementary and unwanted alignments are excluded using SAM flag filtering.

## Author

Mohamed Bahrudeen
Tampere University
