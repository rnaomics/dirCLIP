# Read-wise Mutation Site Co-occurrence Analysis

## Purpose

This script identifies mutation sites that co-occur within individual transcript-aligned reads.

The workflow:

1. Adds MD and NM tags to BAM files.
2. Parses aligned reads.
3. Identifies mismatches and deletions.
4. Extracts mutation-site positions per read.
5. Stores mutation-site coordinates for downstream co-occurrence analyses.

---

## Input Files

### BAM File

```text
filt_q8.bam
```

### Transcript Reference FASTA

```text
gencode.vM1.annotation_transcripts_with_chr.fa
```

---

## Workflow

```text
BAM file
    ↓
Add MD tags
    ↓
Parse aligned reads
    ↓
Detect mismatches and deletions
    ↓
Extract mutation positions
    ↓
Group mutation sites by read
    ↓
Mutation-location table
```

---

## Mutation Identification

The script considers:

- Mismatches identified through lowercase bases in MD-tag-aware alignments.
- Deletions identified through missing read positions in aligned pairs.

For each read, all mutation coordinates are recorded and stored as a comma-separated list.

---

## Parameters

```python
batch_size = 100000
threads = 2
```

---

## Output Files

### MD-tagged BAM

```text
filt_q8_md.bam
```

BAM file containing MD and NM tags.

### Indexed BAM

```text
filt_q8_md.bam.bai
```

Index for the MD-tagged BAM.

### Mutation Locations

```text
filt_q8_mut_loc.tsv.gz
```

Columns:

```text
read_id
ref_name
ref_mod_pos
```

where `ref_mod_pos` contains all mutation positions observed within a read.

---

## Output Structure

```text
Result/5_Co-occurance_sites/
└── 3_SRSF3_DIRCLIP/
    └── Co-SRSF3_dcDNA/
        └── rep1_1pct/
            ├── filt_q8_md.bam
            ├── filt_q8_md.bam.bai
            └── filt_q8_mut_loc.tsv.gz
```

---

## Notes

- Unmapped reads are excluded.
- Mutation coordinates are reported using 1-based reference positions.
- Processing is performed in batches to reduce memory usage.
- The generated mutation-location table is intended for downstream mutation-site co-occurrence analyses.

---

## Author

**Mohamed Bahrudeen**
Tampere University
