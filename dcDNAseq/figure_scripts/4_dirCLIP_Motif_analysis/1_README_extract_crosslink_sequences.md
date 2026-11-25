# Crosslink Site Sequence Extraction

## Purpose
This workflow extracts genomic sequences surrounding dirCLIP crosslinking sites and exports both BED and FASTA files.

## Inputs
- dirCLIP site file: `Pos_mut_thresh_<threshold>.tsv`
  - Column 1: chromosome
  - Column 2: genomic position
  - Column 3: score
- GENCODE mm9 gene annotation GTF
- BSgenome.Mmusculus.UCSC.mm9

## Processing Steps
1. Load dirCLIP sites.
2. Remove sites on random chromosomes.
3. Create a GRanges object from crosslinking positions.
4. Expand each position to a fixed window (`win`).
5. Assign gene strand using overlaps with annotated genes.
6. Remove ambiguous sites overlapping genes on both strands.
7. Export stranded regions as BED.
8. Extract genomic sequences using `getSeq()`.
9. Export sequences as compressed FASTA.

## Window Size
Current setting:

```r
wins <- c(100)
```

A 100-nt sequence centered on each crosslink site is extracted.

## Outputs
### BED
`Pos_mut_thresh_<threshold>.bed.gz`

Contains stranded genomic intervals.

### FASTA
`Pos_mut_thresh_<threshold>.fasta.gz`

Contains extracted sequences suitable for motif analysis tools such as MEME, STREME, or HOMER.

## Running
```bash
Rscript extract_crosslink_sequences.R
```

## Key Packages
- GenomicRanges
- rtracklayer
- Biostrings
- BSgenome.Mmusculus.UCSC.mm9
- dplyr
