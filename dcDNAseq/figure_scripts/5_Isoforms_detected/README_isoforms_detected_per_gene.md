# Isoforms Detected Per Gene Analysis (Read Count and TPM-Based)

## Purpose

This script quantifies the number of detected transcript isoforms per gene and compares isoform diversity between SRSF3 dirCLIP and TotalRNA control samples.

The analysis is performed using both:

- Raw transcript read counts
- TPM-normalized transcript abundances

To focus the analysis on mutation-associated transcript diversity, the workflow:

1. Identifies transcripts containing dirCLIP mutation sites.
2. Resets SRSF3 CLIP transcript counts to zero for transcripts lacking mutation sites.
3. Retains only genes containing at least one mutation-site-associated transcript detected in the SRSF3 CLIP sample.

The resulting analysis therefore focuses on genes associated with identified dirCLIP mutation sites.

---

## Input Files

### Transcript Counts

```text
filt_q8_read_counts_per_transcript.tsv
```

### Significant Mutation Sites

```text
Pos_mut_thresh_0.03.tsv
```

### Transcript Annotation

```text
gencode.vM1.annotation.gtf.gz
```

---

## Workflow

```text
Transcript counts
        ↓
Merge replicates
        ↓
Calculate TPM values
        ↓
Identify transcripts containing dirCLIP mutation sites
        ↓
Reset SRSF3 CLIP counts for transcripts lacking mutation sites
        ↓
Retain genes containing mutation-site-associated transcripts
        ↓
Determine expressed isoforms
        ↓
Count isoforms per gene
        ↓
Generate summary bar plots
```

---

## Analysis

For each gene:

- Expressed transcript isoforms are identified.
- Transcript expression is evaluated using either raw read counts or TPM values.
- The number of expressed isoforms per gene is calculated.
- Distribution of isoform counts across genes is summarized.
- Results are reported using both raw counts and TPM thresholds.

Genes lacking mutation-site-associated transcripts in the SRSF3 CLIP sample are excluded from downstream analyses.
---

## Parameters

```r
Scale_percents = c(100)
Read_count_thresholds = c(5)
tpm_thresh = 0.5
```

---

## Output Figures

### Full Isoform Distribution

```text
N_isoform_per_gene_barplot.pdf
N_isoform_per_gene_barplot_norm.pdf
```

### Grouped Isoform Distribution

```text
more_than_5_isoforms_per_gene_barplot.pdf
more_than_5_isoforms_per_gene_barplot_norm.pdf
```

Grouped categories:

```text
1
2
3
4
≥5
```

---

## Output Directory

```text
Result/figure/8_Isoforms_detected_per_gene/
```

---

## Notes

- TPM values are calculated using transcript lengths derived from GENCODE exon annotations.
- Biological replicates are merged before downstream analyses.
- Transcripts lacking identified dirCLIP mutation sites are assigned zero counts in the SRSF3 CLIP sample before isoform-count calculations.
- Genes are retained only if at least one mutation-site-associated transcript is detected in the SRSF3 CLIP sample.
- Isoform diversity is quantified using both raw read counts and TPM-normalized expression values.

---

## Author

Mohamed Bahrudeen
Tampere University
