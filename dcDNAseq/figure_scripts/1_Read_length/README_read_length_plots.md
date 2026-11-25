# Read Length Visualization

## Purpose

This script generates read-length distributions.

## Input

The script expects `read_length.txt` files generated from BAM alignments.

## Figures Generated

- `RL_violin_from_bam.pdf`
  - Violin plot of read-length distributions.

- `RL_2d_prob_heatmap_from_bam.pdf`
  - 2D probability heatmap showing read-length ranges across samples.

## R Dependencies

```r
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)
```

## Workflow

```text
BAM files
   ↓
read_length.txt
   ↓
Read-length distributions
   ↓
Violin plots
2D probability heatmaps
```

## Output Directory

```text
Results/6_Read_length/
├── RL_violin_from_bam.pdf
└── RL_2d_prob_heatmap_from_bam.pdf
```


## Author

Mohamed Bahrudeen
Tampere University
