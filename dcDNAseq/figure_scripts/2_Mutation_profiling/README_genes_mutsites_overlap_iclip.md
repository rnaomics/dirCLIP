# Genes with Mutation Sites Overlapping iCLIP Peaks

## Purpose

This script quantifies and visualizes the number of genes containing significant mutation sites that overlap iCLIP peaks.

The analysis compares two background conditions:

- TotalRNA cDNA no UV (-UV)
- TotalRNA cDNA UV (+UV)

---

## Input Files

### Gene-level mutation summaries

```text
GW_mut_thresh_0.03.tsv
```

Example:

```text
figure/figure_data/1_Mutation_profiling/
├── ref_TotalRNA_cDNA_no_UV/
└── ref_TotalRNA_cDNA_UV/
```

---

## Workflow

```text
Gene-level mutation summaries
            ↓
Count genes containing mutation sites
            ↓
Determine overlap with iCLIP peaks
            ↓
Classify genes
      ├── Overlap with iCLIP
      └── No overlap with iCLIP
            ↓
Generate stacked bar plot
```

---

## Analysis

For each background condition:

- Genes with at least one mutation site are identified.
- Genes containing mutation sites overlapping iCLIP peaks are counted.
- Remaining genes are classified as non-overlapping.
- Counts are visualized as a stacked bar plot.

---

## Parameters

```r
peak_caller = "iCount"
y_tick_interval = 500
```

---

## Output Figure

```text
SRSF3_Genes_with_Mutsites_overlap_with_iCLIPsites_bar.pdf
```

The figure shows:

- Number of genes with mutation sites overlapping iCLIP peaks.
- Number of genes with mutation sites not overlapping iCLIP peaks.
- Comparison between +UV and -UV controls.

---

## Output Location

```text
Result/figure/1_Mutation_profiling/
└── SRSF3_Genes_with_Mutsites_overlap_with_iCLIPsites_bar.pdf
```

---

## Color Scheme

```r
Colors = c(
  "Yes" = "gray40",
  "No"  = "gray70"
)
```

---

## Notes

- Overlap is determined using the column:

```text
Sites_overlap_with_iCount
```

- Genes are counted once regardless of the number of overlapping sites.
- This is a manuscript-specific figure-generation script and retains fixed sample definitions and plotting parameters.

---

## Author

**Mohamed Bahrudeen**
Tampere University
