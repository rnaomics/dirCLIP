# BAM Merging and Subsampling Workflow

## Purpose

This script merges BAM files from biological replicates and generates subsampled BAM datasets for downstream reproducibility and sensitivity analyses.

The workflow:

1. Merges replicate BAM files.
2. Sorts merged alignments.
3. Generates independent subsampled datasets.
4. Indexes subsampled BAM files.
5. Converts subsampled BAM files to BED format.

---

## Input Files

### BAM files

```text
filt_q8.bam
```

Located in replicate-specific alignment directories.

---

## Workflow

```text
Replicate BAM files
          ↓
Merge BAM files
          ↓
Sort merged BAM
          ↓
Generate subsamples
          ↓
Index BAM files
          ↓
Convert BAM to BED
```

---

## Parameters

```python
threads = 2
N = 3
sample_percent = 70
```

---

## Output Files

### Merged BAM

```text
filt_q8.bam
```

### BAM Index

```text
filt_q8.bam.bai
```

### Subsampled BAM

```text
SS1/filt_q8.bam
SS2/filt_q8.bam
SS3/filt_q8.bam
```

### Subsampled BED

```text
filt_q8.bed.gz
```

Generated from each subsampled BAM file.

---

## Output Structure

```text
Hmga1A_CLIP_cDNA_merge/
└── subsample_70/
    ├── SS1/
    ├── SS2/
    └── SS3/
```

---

## External Tools

- samtools
- bedtools
- pigz

---

## Notes

- Independent random seeds are used for each subsample.
- BAM files are sorted after merging.
- BED files are compressed with pigz.

---

## Author

**Mohamed Bahrudeen**
Tampere University
