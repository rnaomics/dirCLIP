# iCLIP Sequence Extraction and Random Background Generation

## Purpose

This script prepares iCLIP peak sequences and matched random background regions for downstream motif-enrichment analyses.

The workflow:

1. Loads iCLIP peak regions.
2. Extends peaks around the peak center.
3. Annotates genomic regions.
4. Removes intergenic regions.
5. Extracts iCLIP FASTA sequences.
6. Generates GC- and annotation-matched random background regions.
7. Exports BED and FASTA files for motif analysis.

---

## Input Files

### iCLIP Peaks

```text
figure/figure_data/2_iclip/SRSF3/icount/
└── SRSF3_merged.3nt_3nt.peaks.bed.gz
```

### Genome Annotation

```text
genome/mm9/Annotation/gencode.vM1.annotation.gtf
```

### Reference Genome

```text
BSgenome.Mmusculus.UCSC.mm9
```

---

## Workflow

```text
iCLIP peaks
      ↓
Peak-center extension
      ↓
Region annotation
      ↓
Remove intergenic peaks
      ↓
Extract genomic sequences
      ↓
Generate GC-matched controls
      ↓
Export BED and FASTA files
```

---

## Parameters

```r
win = 100
```

### Peak Centering

```r
From_peak_centre = TRUE
```

### Random Background Replicates

```r
N = 10
```

### GC-Matching Draws

```r
n_draws = 100
```

### Peak Extension Window

```r
Left_offset  = 50
Right_offset = 50
```

Peaks are expanded to 100-nt windows centered on the peak midpoint.

---

## Output Structure

```text
Result/figure/7_iclip/
└── SRSF3/
    └── icount/
        └── Motif_analysis/
            └── win_100/
                ├── iCLIP.bed.gz
                ├── iCLIP.fasta.gz
                └── random_region/
                    ├── random_region_1.bed.gz
                    ├── random_region_1.fasta.gz
                    ├── random_region_2.bed.gz
                    ├── random_region_2.fasta.gz
                    └── ...
```

---

## Output Files

### iCLIP Peak Regions

```text
iCLIP.bed.gz
```

Extended iCLIP peak regions after filtering and annotation.

### iCLIP FASTA Sequences

```text
iCLIP.fasta.gz
```

Sequences extracted from the extended iCLIP peak regions.

### Random Background BED Files

```text
random_region_*.bed.gz
```

GC- and annotation-matched random genomic regions.

### Random Background FASTA Files

```text
random_region_*.fasta.gz
```

Sequences extracted*from the matched random background*regions.

---

## Region Matching Strategy

Random background regions are matched using:

- Genomic annotation class
  - 5'UTR
  - CDS
  - 3'UTR
  - Exon
  - Intron

- GC content

Intergenic regions are excluded.

---

## Notes

- Chromosomes containing the term `random` are removed.
- Peak regions are centered before extension.
- GC matching is performed independently for each genomic feature class.
- The generated random regions are intended for downstream motif-enrichment analysis.

---

## Author

**Mohamed Bahrudeen**
Postdoctoral Research Fellow
Tampere University
