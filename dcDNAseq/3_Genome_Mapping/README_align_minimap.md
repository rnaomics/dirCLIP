# 🧬 Minimap2 Alignment and BAM Processing Pipeline

This pipeline performs alignment and post-processing of Oxford Nanopore sequencing data using minimap2, samtools, and bedtools.

The workflow is modular and uses subcommands, allowing you to run alignment and filtering steps independently.

---

## 📦 Requirements

- Python ≥ 3.7
- minimap2
- samtools
- bedtools
- pysam

Install pysam (if not already installed):

```bash
pip install pysam
```

---

## 📂 Input

Input data should be FASTQ or FASTQ.GZ files:

You will provide FASTQ file paths directly via command-line arguments.
---

## 🚀 Usage

    python align_minimap.py <command> [options]

Available commands:
- align → perform alignment and BAM processing
- filter → filter BAM files by read quality

---

## 1. Alignment

```bash
python align_minimap.py align \
  --fastq_dirs sample1.fastq.gz sample2.fastq.gz \
  --output_dirs out1 out2 \
  --genome genome.fa \
  --splice_sites splice_sites.bed \
  --threads 6
```

---

## 2. Filtering

```bash
python align_minimap.py filter \
  --output_dirs out1 out2 \
  --min_quality 8
```

---

## ⚙️ Parameters

### Align
- --fastq_dirs : input FASTQ files (required)
- --output_dirs : output directories (required)
- --genome : reference genome (required)
- --splice_sites : Splice junction BED file (optional)
- --bam_name : output BAM name (default: sorted.bam)
- --threads : number of threads (default: 4)

### Filter
- --output_dirs : directories containing BAM files
- --bam_name : input BAM name (default: sorted.bam)
- --min_quality : minimum mean read quality (default: 0)

---
## ▶️ Example Run (Step-by-step)

### 1. Go to working directory

```bash
cd Tools/dcDNAseq
```

### 2. Run alignment

```bash
python 3_Genome_Mapping/align_minimap.py align \
  --fastq_dirs \
  ../1_Test_data/SRSF3_dcDNA/rep1_1pct.fastq.gz \
  ../1_Test_data/SRSF3_dcDNA/rep2_1pct.fastq.gz \
  ../1_Test_data/SRSF3_dcDNA/rep3_1pct.fastq.gz \
  ../1_Test_data/TotalRNA_dcDNA_no_UV/rep1_0.2pct.fastq.gz \
  ../1_Test_data/TotalRNA_dcDNA_no_UV/rep2_0.2pct.fastq.gz \
  ../1_Test_data/TotalRNA_dcDNA_no_UV/rep3_0.2pct.fastq.gz \
  --output_dirs \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
  --genome ../genome/mm9/mm9.fa \
  --splice_sites ../genome/mm9/Annotation/mm9.knownGene.ss.bed \
  --threads 2
```

## Optional: Run without splice sites

```bash
python 3_Genome_Mapping/align_minimap.py align \
  --fastq_dirs \
  ../1_Test_data/SRSF3_dcDNA/rep1_1pct.fastq.gz \
  ../1_Test_data/SRSF3_dcDNA/rep2_1pct.fastq.gz \
  ../1_Test_data/SRSF3_dcDNA/rep3_1pct.fastq.gz \
  ../1_Test_data/TotalRNA_dcDNA_no_UV/rep1_0.2pct.fastq.gz \
  ../1_Test_data/TotalRNA_dcDNA_no_UV/rep2_0.2pct.fastq.gz \
  ../1_Test_data/TotalRNA_dcDNA_no_UV/rep3_0.2pct.fastq.gz \
  --output_dirs \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
  --genome ../genome/mm9/mm9.fa \
  --threads 2
```bash

### 3. Run filtering
```bash
python 3_Genome_Mapping/align_minimap.py filter \
  --output_dirs \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/ \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/ \
  ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
  ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
  --min_quality 8
```
---
## 📤 Output

Each output directory contains:

```
output/
├── sorted.bam
├── sorted.bam.bai
├── sorted_alignmentStats.txt
├── sorted_flag_stats.txt
├── sorted.bed.gz
├── filt_qX.bam        (if filtering enabled)
├── filt_qX.bam.bai
├── filt_qX_alignmentStats.txt
├── filt_qX.bed.gz
```

---

## 🔁 Workflow

Align:
1. Run minimap2
2. Convert and sort BAM
3. Index BAM
4. Generate stats
5. Convert to BED

Filter:
1. Compute mean read quality
2. Filter reads
3. Write filtered BAM
4. Generate stats and BED

---

## 📌 Notes

- Input and output directories must be paired
- Filtering is skipped if `--min_quality 0`
- Required tools must be available in system PATH

---
## 👨‍🔬 Author

Mohamed Mohamed Bahrudeen
