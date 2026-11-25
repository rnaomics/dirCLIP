# Mutation Sites Per Read and Inter-site Distance Analysis

## Purpose

This script quantifies significant dirCLIP mutation sites at the read level and isoform level, and characterizes the spacing between co-occurring mutation sites within individual reads.

The workflow produces:

1. Mutation-site counts per read.
2. Mutation-site counts per isoform.
3. Distance distributions between mutation sites within reads.

---

## Input Files

### Read-level mutation locations

```text
filt_q8_mut_loc.tsv.gz
```

Generated from the mutation-site co-occurrence workflow.

### Significant mutation sites

```text
Pos_mut_thresh_0.03.tsv
```

---

## Workflow

```text
Read-level mutation locations
              ↓
Significant mutation-site catalog
              ↓
Assign mutation sites to reads
              ↓
Count mutation sites per read
              ↓
Count mutation sites per isoform
              ↓
Calculate distances between mutation sites
              ↓
Generate summary tables and figures
```

---

## Analysis

### Mutation Sites Per Read

For each read:

- Mutation positions are matched against the significant mutation-site catalog.
- The number of significant mutation sites within the read is counted.
- Reads are grouped according to the number of mutation sites.

### Mutation Sites Per Isoform

For each transcript:

- The total number of significant mutation sites is calculated.

### Distances Between Mutation Sites

For reads containing multiple mutation sites, the script calculates:

- Distance from the first mutation site.
- Distance from the last mutation site.
- Distance between consecutive mutation sites.

Both raw counts and weighted counts are computed.

---

## Distance Bins

```text
1-25 bp
26-100 bp
101-200 bp
201-500 bp
>500 bp
```

---

## Output Files

### Mutation-site summary table

```text
df_sites_per_read.tsv.gz
```

Contains:

```text
read_id
ref_name
ref_mut_pos
N_mut_per_read
N_mut
```

### Mutation Sites per Read Figure

```text
mut_sites_per_read.pdf
```

Contains:

- Distribution of mutation sites per read.
- Distribution of mutation sites per isoform.

### Mutation-site Distance Figure

```text
dist_mut_sites_per_read.pdf
```

Contains:

- Distance from first mutation site.
- Distance from last mutation site.
- Distance between consecutive mutation sites.
- Raw-count and weighted-count summaries.

---

## Output Directory

```text
Result/figure/
└── 5_Co-occurance_sites/
    └── 3_SRSF3_DIRCLIP/
        └── SRSF3_dcDNA/
            └── Figures/
                ├── mut_sites_per_read.pdf
                └── dist_mut_sites_per_read.pdf
```

---

## Notes

- Only mutation positions matching the significant mutation-site catalog are retained.
- Reads aligning to transcripts without significant mutation sites are excluded.
- Weighted counts normalize for the number of distance observations contributed by each read.

---

## Author

**Mohamed Bahrudeen**  
Tampere University
