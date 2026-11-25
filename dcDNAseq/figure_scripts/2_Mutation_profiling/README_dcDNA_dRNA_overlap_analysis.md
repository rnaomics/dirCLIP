# dcDNA and dRNA dirCLIP Overlap Analysis

## Purpose

This script compares dirCLIP sites or regions identified from dcDNA-based and dRNA-based analyses.

The workflow:

1. Converts transcript coordinates to genomic coordinates.
2. Maps dcDNA and dRNA dirCLIP signals to a common genome assembly.
3. Performs overlap analysis between datasets.
4. Computes overlap statistics.
5. Generates UpSet plots summarizing shared and unique sites.

---

## Input Files

### dcDNA dirCLIP Sites or Regions

```text
Pos_mut_thresh_0.03.tsv
Pos_mut_thresh_0.03_region.tsv
```

### dRNA dirCLIP Results

```text
Sup_table1.xlsx
Sup_table3.xlsx
```

### Genome Annotations

```text
gencode.vM1.annotation.gtf.gz
gencode.vM36.annotation.gtf.gz
```

### Liftover Chain

```text
mm39ToMm9.over.chain
```

---

## Workflow

```text
dcDNA dirCLIP
        ↓
Transcript-to-genome conversion
        ↓
Liftover between assemblies
        ↓
dRNA dirCLIP genomic sites
        ↓
Overlap analysis
        ↓
Statistics table
        ↓
UpSet plot
```

---

## Analysis Outputs

### Overlap Statistics

Reported metrics include:

- Total dcDNA sites
- Total dRNA sites
- Overlapping dcDNA sites
- Overlapping dRNA sites
- Percent overlap for each dataset

### UpSet Plot

Visualizes:

- Shared sites/regions
- dcDNA-specific sites/regions
- dRNA-specific sites/regions

---

## Output Files

```text
*_olap_with_dRNA_sites.tsv
*_olap_with_dRNA_regions.tsv
```

```text
dcDNA_dRNA_sites_overlap.pdf
dcDNA_dRNA_regions_overlap.pdf
```

---

## Notes

- Transcript coordinates are converted using TxDb transcript models.
- Liftover is used to compare datasets generated on different mouse genome assemblies.
- Site widths are expanded to 11 bp before overlap analysis.
- This is a manuscript-specific comparison workflow and retains fixed parameters and thresholds.

---

## Author

Mohamed Bahrudeen
Tampere University
