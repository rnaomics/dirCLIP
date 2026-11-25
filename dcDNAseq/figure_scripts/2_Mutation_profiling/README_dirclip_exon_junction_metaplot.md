# Distance of dirCLIP Mutation Sites from Exon Junctions

## Purpose

This script analyzes the positional distribution of significant dirCLIP mutation sites relative to exon donor and acceptor splice sites.

The workflow determines the nearest exon junction for each mutation site and generates density plots showing enrichment around splice junctions.

---

## Input Files

### Significant Mutation Sites

```text
Pos_mut_thresh_0.03.tsv
```

### Gene Annotation

```text
genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz
```

---

## Workflow

```text
Significant mutation sites
            ↓
Load exon annotations
            ↓
Identify donor sites
            ↓
Identify acceptor sites
            ↓
Calculate strand-aware distances
            ↓
Generate density plots
            ↓
Metaplot around exon junctions
```

---

## Analysis

For each mutation site:

- Nearest exon donor site is identified.
- Nearest exon acceptor site is identified.
- Strand-aware distances are calculated.
- Density distributions are generated for donor and acceptor sites separately.

---

## Parameters

```r
dist_thresh = 200
```

Only sites within ±200 bp are shown in the final metaplots.

---

## Output Figure

```text
SRSF3_mutsites_distance_from_exon_ends_metaplot.pdf
```

The figure contains:

- Distance distribution from exon donor sites.
- Distance distribution from exon acceptor sites.

---

## Output Location

```text
Result/figure/1_Mutation_profiling/
└── SRSF3_mutsites_distance_from_exon_ends_metaplot.pdf
```

---

## Notes

- Distances are calculated in a strand-aware manner.
- Donor and acceptor junctions are defined according to transcript strand.
- Exon annotations are derived from Gencode.
- Density plots are generated using ggplot2.

---

## Author

**Mohamed Bahrudeen**
Tampere University
