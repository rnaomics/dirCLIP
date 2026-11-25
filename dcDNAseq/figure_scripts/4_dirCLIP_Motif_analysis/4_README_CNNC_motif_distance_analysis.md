# CNNC Distance Analysis Around Mutation Sites

## Purpose

This script quantifies the distance of mutation sites to nearby CNNC motifs and compares the observed distribution against two control strategies:

1. Sequence shuffling controls.
2. GC-matched random-background sequence controls.

The resulting distance distributions can be used to evaluate whether CNNC motifs preferentially occur near mutation sites.

---

## Input Files

### Mutation-site FASTA sequences

```text
Pos_mut_thresh_0.03.fasta.gz
```

### Random background FASTA sequences

```text
rnd_bg_1.fasta.gz
rnd_bg_2.fasta.gz
rnd_bg_3.fasta.gz
```

---

## Workflow

```text
Mutation-site FASTA
          ↓
Extract k-mers
          ↓
Identify CNNC motifs
          ↓
Compute nearest CNNC distance
          ↓
Generate distance distribution
          ↓
Compare with:
   • Shuffled sequences
   • Random-background sequences
```

---

## CNNC Motifs

```r
CNNC_motifs = c(
  "TCTTC", "CTTCT", "CTTCC", "CTTCA", "CAGCA",
  "CTCCT", "CATCT", "TCCTC", "CCTTC", "CCTCC"
)
```

---

## Parameters

```r
Coverage_thresholds = c(50)
Diff_thresholds = c(0.03)
kmers = c(5)
window_length = 100
N = 3
```

---

## Output Files

### Observed Distribution

```text
CNNC_dist_thresh_0.03.tsv
```

### Shuffled Controls

```text
CNNC_dist_randomers_thresh_0.03.tsv
```

### Random Background Controls

```text
CNNC_dist_rnd_bg_thresh_0.03.tsv
```

Each table contains CNNC distances and corresponding frequencies.

---

## Notes

- Distances are measured relative to the center of the sequence window.
- Sequence shuffling preserves nucleotide composition while disrupting motif organization.
- Random-background controls use independently generated background sequences.
- This is a manuscript-specific support script and retains fixed paths and parameters.

---

## Author

**Mohamed Bahrudeen**
Tampere University
