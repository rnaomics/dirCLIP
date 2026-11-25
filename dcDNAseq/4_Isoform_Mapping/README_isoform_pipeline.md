# Isoform-Specific Read Assignment and Alignment Pipeline

This pipeline performs isoform-specific read assignment and transcript-level alignment of Oxford Nanopore sequencing reads.

## Workflow

1. Transcriptome preparation from GTF annotation
2. Bambu-based read-to-isoform assignment
3. FASTQ sorting by isoform
4. Isoform-specific alignment using minimap2
5. Transcript-level read counting

---

# Requirements

## Software

- Python ≥ 3.8
- R
- Bambu
- minimap2
- samtools
- bedtools
- gffread
- pigz

## Python Packages

```bash
pip install pandas numpy pysam biopython
```

## Script path

```bash
cd Tool_submission/4_Isoform_Mapping
```

---

# 1. Prepare Transcriptome

```bash
python isoform_pipeline.py prepare \
--gtf ../genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz \
--genome_fasta ../genome/mm9/mm9.fa \
--transcript_dir ../genome/mm9/Annotation/Transcripts/
```

# 2. Read Assignment

```bash
python isoform_pipeline.py assign \
--bam_files \
    ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/filt_q8.bam \
--output_dirs \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--genome_fasta ../genome/mm9/mm9.fa \
--gtf ../genome/mm9/Annotation/gencode.vM1.annotation_with_chr.gtf.gz
```

# 3. Sort FASTQ by Isoform

```bash
python isoform_pipeline.py sort-fastq \
--fastq_files \
    ../1_Test_data/SRSF3_dcDNA/rep1_1pct.fastq.gz \
    ../1_Test_data/SRSF3_dcDNA/rep2_1pct.fastq.gz \
    ../1_Test_data/SRSF3_dcDNA/rep3_1pct.fastq.gz \
    ../1_Test_data/TotalRNA_dcDNA_no_UV/rep1_0.2pct.fastq.gz \
    ../1_Test_data/TotalRNA_dcDNA_no_UV/rep2_0.2pct.fastq.gz \
    ../1_Test_data/TotalRNA_dcDNA_no_UV/rep3_0.2pct.fastq.gz \
--output_dirs \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
```

# 4. Isoform-Specific Alignment

```bash
python isoform_pipeline.py align \
--output_dirs \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--transcript_dir ../genome/mm9/Annotation/Transcripts/ \
--bam_name filt_q8.bam \
--threads 2 \
--sam_flag 2308 \
--delete_intermediate
```

# 5. Count Reads per Transcript/Gene

```bash
python isoform_pipeline.py count \
--output_dirs \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
--bam_name filt_q8.bam \
--gtf ../genome/mm9/Annotation/gencode.vM1.annotation_with_chr.gtf.gz \
--countby gene_name
```

# 6. Run Complete Workflow

```bash
python isoform_pipeline.py all \
    --bam_files \
    ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/filt_q8.bam \
    ../Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/filt_q8.bam \
    --fastq_files \
    ../1_Test_data/SRSF3_dcDNA/rep1_1pct.fastq.gz \
    ../1_Test_data/SRSF3_dcDNA/rep2_1pct.fastq.gz \
    ../1_Test_data/SRSF3_dcDNA/rep3_1pct.fastq.gz \
    ../1_Test_data/TotalRNA_dcDNA_no_UV/rep1_0.2pct.fastq.gz \
    ../1_Test_data/TotalRNA_dcDNA_no_UV/rep2_0.2pct.fastq.gz \
    ../1_Test_data/TotalRNA_dcDNA_no_UV/rep3_0.2pct.fastq.gz \
    --output_dirs \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep1_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep2_1pct/ \
    ../Result/4_Isoform_Mapping/SRSF3_dcDNA/rep3_1pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/ \
    ../Result/4_Isoform_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/ \
    --genome_fasta ../genome/mm9/mm9.fa \
    --gtf ../genome/mm9/Annotation/gencode.vM1.annotation_with_chr.gtf.gz \
    --transcript_dir ../genome/mm9/Annotation/Transcripts/ \
    --bam_name filt_q8.bam \
    --threads 2 \
    --sam_flag 2308 \
    --countby gene_name \
    --delete_intermediate
```

---

# Output Files

```text
sample/
├── Out/
│   └── readToTranscriptMaps.tsv
├── isoform_counts.csv
├── filt_q8.bam
├── filt_q8.bam.bai
├── filt_q8.bed.gz
├── filt_q8_alignmentStats.txt
├── filt_q8_read_counts_per_transcript.tsv
└── fastq_by_isoform/
```

---

# Notes

- Transcript IDs are prefixed with `Chr_` during transcriptome preparation.
- Only uniquely assigned reads are used for isoform-specific FASTQ generation.
- `--delete_intermediate` removes temporary files after successful completion.
- Input BAM files should be sorted and indexed.

---

