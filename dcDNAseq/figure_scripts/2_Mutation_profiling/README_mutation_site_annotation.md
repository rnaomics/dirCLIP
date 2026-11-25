# Mutation Site Annotation Analysis

## Purpose

This script annotates significant mutation sites with genomic features and visualizes their distribution across transcript annotations.

The analysis classifies mutation sites into:

- 3UTR
- 5UTR
- Exon
- Intron
- Intergenic

and summarizes the fraction of sites in each category.

---

## Input Files

### Significant Mutation Sites

```text
Pos_mut_thresh_0.03.tsv
```

### Genome Annotation

```text
gencode.vM1.annotation.gtf.gz
```

---

## Workflow

```text
Mutation sites
      ↓
Load GTF annotation
      ↓
Generate exon regions
      ↓
Generate intron regions
      ↓
Identify 5'UTRs and 3'UTRs
      ↓
Annotate sites
      ↓
Calculate annotation fractions
      ↓
Generate stacked bar plot
```

---

## Annotation Categories

```text
3UTR
5UTR
Exon
Intron
Intergenic
```

Gene-type annotations are also assigned when available.

---

## Output Figure

```text
Hmga1a_1b_Mut_sites_annotation_stacked_bar.pdf
```

The figure displays the fraction of mutation sites overlapping each genomic annotation.

---

## Output Location

```text
Result/figure/1_Mutation_profiling/
└── Hmga1a_1b_Mut_sites_annotation_stacked_bar.pdf
```

---

## Intermediate Annotation Files

```text
gencode_exons_and_intron.RData
gencode_UTR.RData
```

---

## Notes

- Random chromosomes are excluded.
- Introns are derived using transcript minus exon intervals.
- UTRs are classified using start and stop codons.
- Site annotations are assigned by genomic overlap.
- This is a manuscript-specific figure-generation script and retains fixed parameters and paths.

---

## Author

**Mohamed Bahrudeen**
Postdoctoral Research Fellow
Tampere University
