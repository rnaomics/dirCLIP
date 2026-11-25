# Mutation Site Distribution Across Transcript Length

## Purpose

This script visualizes the distribution of significant dirCLIP mutation sites across transcript bodies after normalizing coordinates by transcript length.

The analysis evaluates whether mutation sites preferentially occur at specific relative positions along transcripts.

---

## Input Files

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
Mutation sites
        ↓
Load transcript annotation
        ↓
Calculate transcript lengths
        ↓
Normalize mutation positions
        ↓
Density estimation
        ↓
Generate distribution plots
```

---

## Analysis

For each mutation site:

- Transcript coordinates are extracted.
- Transcript length is calculated from exon annotations.
- Mutation position is normalized:

```text
normalized_position = position / transcript_length
```

Values range from:

```text
0 = transcript start
1 = transcript end
```

---

## Output Figure

```text
SRSF3_Mutsites_distribution_across_length.pdf
```

Contains:

- Filled density distribution.
- Line density distribution.

Both panels display mutation-site density across normalized transcript coordinates.

---

## Output Location

```text
Result/figure/1_Mutation_profiling/
└── SRSF3_Mutsites_distribution_across_length.pdf
```

---

## Notes

- Transcript lengths are calculated by summing exon lengths.
- Transcript identifiers are matched to mutation-site annotations.
- Coordinates are normalized to enable comparison across transcripts of different lengths.

---

## Author

**Mohamed Bahrudeen**
Postdoctoral Research Fellow
Tampere University
