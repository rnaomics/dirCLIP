# dcDNA-seq analysis Folder

## Overview

This folder contains the computational workflows, analysis scripts,
example sequencing data, processed analysis tables, and
figure-generation scripts associated with the manuscript:

**dirCLIP PROFILES VARIANT-SPECIFIC RNA-PROTEIN INTERACTIONS VIA NANOPORE LONG-READ SEQUENCING**

The code supports processing Oxford Nanopore sequencing data and
downstream analyses of genome- and transcript/isoform-aligned reads,
mutation profiles, crosslink-associated sequence motifs, isoform
detection, and mutation-site co-occurrence. The `figure_scripts/`
directory contains scripts used for manuscript figures and related
analyses.

The folder is organised so that readers can distinguish the main
data-processing workflows from the scripts used to generate manuscript
figures. Detailed, script-specific instructions are provided in README
files within the relevant directories.

> **Important:** The folder includes small, subsampled FASTQ example
> datasets and processed figure data. It does not include the full raw
> sequencing runs or all external reference resources used by the
> workflows. Full reproduction therefore requires access to the
> appropriate input data, genome/transcriptome reference files,
> annotations, and software listed in the relevant workflow
> documentation.

## Folder contents

``` text
dcDNAseq/
├── 1_Test_data/
│   ├── SRSF3_dcDNA/
│   └── TotalRNA_dcDNA_no_UV/
├── 2_Basecalling/
│   ├── dorado_basecall/
│   └── guppy_basecall/
├── 3_Genome_Mapping/
├── 4_Isoform_Mapping/
├── 5_Mutation_profiling/
├── figure_scripts/
│   ├── 1_Read_length/
│   ├── 2_Mutation_profiling/
│   ├── 3_iCLIP_Motif_analysis/
│   ├── 4_dirCLIP_Motif_analysis/
│   ├── 5_Isoforms_detected/
│   ├── 6_Co-occurence_analysis/
│   └── figure_data/
└── README.md
```

The directory names are retained as they appear in the supplied
folder.

## 1. Test data

**Directory:** `1_Test_data/`

Contains small, subsampled FASTQ datasets intended to let users inspect
the input format and test parts of the workflows without downloading or
processing the complete sequencing datasets.

  -----------------------------------------------------------------------
  Subdirectory                        Contents
  ----------------------------------- -----------------------------------
  `SRSF3_dcDNA/`                      Three example replicate FASTQ files
                                      (`rep1_1pct.fastq.gz`,
                                      `rep2_1pct.fastq.gz`, and
                                      `rep3_1pct.fastq.gz`) for the SRSF3
                                      dcDNA sample.

  `TotalRNA_dcDNA_no_UV/`             Three example replicate FASTQ files
                                      (`rep1_0.2pct.fastq.gz`,
                                      `rep2_0.2pct.fastq.gz`, and
                                      `rep3_0.2pct.fastq.gz`) for the
                                      TotalRNA dcDNA no-UV control.
  -----------------------------------------------------------------------

The percentage in each filename indicates the subsampling level
represented by the example file. These files are test inputs, not a
replacement for the complete sequencing datasets used in the manuscript.

## 2. Basecalling

**Directory:** `2_Basecalling/`

Contains Python scripts for basecalling Oxford Nanopore sequencing data.

  ----------------------------------------------------------------------------
  File or directory                        Purpose
  ---------------------------------------- -----------------------------------
  `dorado_basecall/dorado_basecaller.py`   Runs basecalling using ONT Dorado.

  `dorado_basecall/README.md`              Dorado-specific usage instructions.

  `guppy_basecall/guppy_basecall.py`       Runs basecalling using ONT Guppy.

  `guppy_basecall/README.md`               Guppy-specific usage instructions.
  ----------------------------------------------------------------------------

**Input:** ONT POD5 or FAST5 data, as supported by the selected
basecaller and its configuration.

**Output:** Basecalled FASTQ reads.

See the README in the corresponding basecaller directory for
command-line options and usage details.

## 3. Genome mapping

**Directory:** `3_Genome_Mapping/`

Contains scripts for aligning basecalled reads to the reference genome
and preparing BAM files for downstream analyses.

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `align_minimap.py`                  Aligns Nanopore reads using
                                      minimap2 and carries out the
                                      alignment-processing steps
                                      described in its accompanying
                                      documentation.

  `bam_merge_and_subsample.py`        Merges replicate BAM files, sorts
                                      and indexes alignments, generates
                                      subsampled BAM datasets, and
                                      converts subsampled alignments to
                                      BED format.
  -----------------------------------------------------------------------

Additional documentation:

-   [`README_align_minimap.md`](3_Genome_Mapping/README_align_minimap.md)
-   [`README_bam_merge_and_subsample.md`](3_Genome_Mapping/README_bam_merge_and_subsample.md)

**Main input:** Basecalled FASTQ files.

**Main output:** Genome-aligned, processed BAM files and, where
requested, subsampled BAM/BED files.

**Main external tools:** minimap2, samtools, bedtools, and pigz (for the
documented BAM-processing workflow).

## 4. Isoform mapping

**Directory:** `4_Isoform_Mapping/`

Contains the workflow for assigning genome-aligned reads to transcript
isoforms and performing transcript/isoform-specific alignment and read
counting.

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `bambu_for_read_assignment.R`       R wrapper for Bambu-based
                                      read-to-transcript assignment and
                                      export of read-to-transcript
                                      mappings.

  `isoform_pipeline.py`               Coordinates transcriptome
                                      preparation, read assignment, FASTQ
                                      sorting by isoform,
                                      isoform-specific alignment, and
                                      transcript-level read counting.
  -----------------------------------------------------------------------

The Python pipeline provides the following workflow stages:

1.  `prepare` -- prepares transcript sequences from the reference genome
    and GTF annotation.
2.  `assign` -- performs Bambu-based read assignment and exports
    read-to-transcript mappings.
3.  `sort-fastq` -- selects uniquely assigned reads and separates FASTQ
    reads by transcript isoform.
4.  `align` -- aligns isoform-specific reads using minimap2 and
    processes the resulting BAM files.
5.  `count` -- calculates read counts per transcript or gene.
6.  `all` -- runs the complete workflow.

Typical outputs documented for this workflow include
`readToTranscriptMaps.tsv`, `isoform_counts.csv`, processed BAM files,
and `filt_q8_read_counts_per_transcript.tsv`.

**Main external tools/packages:** R, Bambu, minimap2, samtools,
bedtools, gffread, and pigz, in addition to the Python packages listed
in the workflow documentation.

See
[`README_isoform_pipeline.md`](4_Isoform_Mapping/README_isoform_pipeline.md)
for the command-line examples, required inputs, parameters, and output
structure.

## 5. Mutation profiling

**Directory:** `5_Mutation_profiling/`

Contains the scripts for per-base mutation profiling and identification
of significant mutation sites.

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `mutation_profiling_pipeline.py`    Runs per-base profiling and
                                      calculates mutation rates from BAM
                                      files. Supports genome-based and
                                      isoform-specific input workflows,
                                      and analyses with or without
                                      deletion inclusion.

  `significant_mutations.R`           Compares mutation rates between
                                      target and control replicates,
                                      applies coverage and
                                      mutation-difference thresholds,
                                      performs statistical testing and
                                      Benjamini--Hochberg FDR correction,
                                      and exports significant sites and
                                      grouped regions.
  -----------------------------------------------------------------------

The mutation-profiling workflow uses the `perbase` tool to obtain
per-position nucleotide counts and calculate mutation rates. The
significant-mutation workflow uses replicate-level mutation-rate files
to identify sites with elevated mutation rates in the target relative to
the control.

The significant-mutation script can generate position-level results, WIG
tracks, and region-level summaries. The documented outputs include files
such as:

-   `Pos_mut_thresh_*.tsv` -- filtered significant mutation sites.
-   `Pos_mut_thresh_*.wig` -- genome-browser tracks, when requested.
-   `Pos_mut_thresh_*_region.tsv` -- grouped mutation-enriched regions.
-   `significant_mutsites.tsv` -- statistical summary of tested sites.

The workflows distinguish analyses that exclude deletions from those
that include deletions. This distinction is reflected in the output
directory structure (`1_No_del` and `2_Inc_del`).

See:

-   [`README_Mutation_Profiling.md`](5_Mutation_profiling/README_Mutation_Profiling.md)
-   [`README_Significant_Mutations.md`](5_Mutation_profiling/README_Significant_Mutations.md)

## 6. Manuscript figure and downstream analysis scripts

**Directory:** `figure_scripts/`

Contains the R and Python scripts used for manuscript figure generation
and related downstream analyses. The processed inputs used by these
scripts are organised in `figure_scripts/figure_data/`.

### 6.1 Read-length analysis

**Directory:** `figure_scripts/1_Read_length/`

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `read_length_statistics.py`         Calculates read-length summary
                                      statistics.

  `read_length_plots.R`               Generates read-length distribution
                                      plots.
  -----------------------------------------------------------------------

See `README_read_length_statistics.md` and `README_read_length_plots.md`
for input and output details.

### 6.2 Mutation-profiling figures

**Directory:** `figure_scripts/2_Mutation_profiling/`

  ------------------------------------------------------------------------------------------
  Script                                                 Purpose
  ------------------------------------------------------ -----------------------------------
  `1_mutation_rate_boxplots.R`                           Plots mutation-rate distributions
                                                         and compares analyses with and
                                                         without deletion inclusion.

  `2_mutation_overlap_upset_plots.R`                     Generates UpSet plots showing
                                                         overlaps among significant mutation
                                                         sites and mutation-enriched
                                                         regions.

  `3_genewise_mutation_sites.R`                          Summarises significant mutation
                                                         sites and regions at the gene
                                                         level, including overlap with iCLIP
                                                         peaks.

  `4_genes_mutsites_overlap_iclip.R`                     Compares genes containing mutation
                                                         sites that overlap iCLIP peaks with
                                                         genes without such overlap.

  `5_dirclip_exon_junction_metaplot.R`                   Analyses the distance of mutation
                                                         sites from exon donor and acceptor
                                                         junctions and generates metaplots.

  `6_mutsites_distribution_across_transcript_length.R`   Plots mutation-site distributions
                                                         across transcript bodies after
                                                         normalising positions by transcript
                                                         length.

  `7_mutation_site_annotation.R`                         Annotates mutation sites using
                                                         genomic features, including UTR,
                                                         exon, intron, and intergenic
                                                         regions.

  `8_dcDNA_dRNA_overlap_analysis.R`                      Performs the dcDNA/dRNA overlap
                                                         analysis described in its
                                                         accompanying README.
  ------------------------------------------------------------------------------------------

Individual README files in this directory describe the corresponding
input tables, analysis steps, parameters, and figure outputs.

### 6.3 iCLIP motif analysis

**Directory:** `figure_scripts/3_iCLIP_Motif_analysis/`

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `extract_iclip_sequence.R`          Extracts sequences around iCLIP
                                      sites for downstream motif
                                      analysis.

  `motif_enrichment_analysis.R`       Performs motif-enrichment analysis
                                      on the extracted sequences.
  -----------------------------------------------------------------------

See `1_README_extract_iclip_sequence.md` and
`2_README_motif_enrichment_analysis.md`.

### 6.4 dirCLIP motif analysis

**Directory:** `figure_scripts/4_dirCLIP_Motif_analysis/`

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `extract_crosslink_sequences.R`     Extracts sequences around dirCLIP
                                      crosslink sites.

  `extract_random_background.R`       Extracts random background
                                      sequences for motif comparisons.

  `motif_enrichment_analysis.R`       Performs motif-enrichment analysis.

  `cnnc_motif_distance_analysis.R`    Analyses the distance of CNNC
                                      motifs relative to crosslink sites.

  `plot_CNNC_distance.R`              Generates plots of the CNNC
                                      motif-distance analysis.
  -----------------------------------------------------------------------

See the numbered README files in this directory for the input
definitions, background extraction, motif analysis, and CNNC distance
workflow.

### 6.5 Isoforms detected per gene

**Directory:** `figure_scripts/5_Isoforms_detected/`

  -----------------------------------------------------------------------
  Script                              Purpose
  ----------------------------------- -----------------------------------
  `isoforms_detected_per_gene.R`      Quantifies detected transcript
                                      isoforms per gene and compares
                                      isoform detection using read-count
                                      and TPM-based approaches.

  -----------------------------------------------------------------------

The documented analysis uses transcript-count tables, mutation-site
information, and transcript annotations. See
`README_isoforms_detected_per_gene.md`.

### 6.6 Mutation-site co-occurrence analysis

**Directory:** `figure_scripts/6_Co-occurence_analysis/`

  ------------------------------------------------------------------------
  Script                               Purpose
  ------------------------------------ -----------------------------------
  `1_cooccurrence_sites.py`            Processes BAM alignments,
                                       identifies mismatch and deletion
                                       positions within reads, and creates
                                       read-level mutation-location
                                       tables.

  `2_plot_mutation_sites_per_read.R`   Summarises significant mutation
                                       sites per read and per isoform,
                                       calculates inter-site distances,
                                       and generates plots.
  ------------------------------------------------------------------------

See `README_cooccurrence_sites.md` and
`README_plot_mutation_sites_per_read.md`.

## 7. Processed data for figure generation

**Directory:** `figure_scripts/figure_data/`

Contains processed tables and supporting datasets used by the manuscript
figure-generation scripts. These are analysis inputs and
intermediate/processed outputs, rather than the primary raw sequencing
data.

The directory is organised by analysis type. The supplied archive
includes, among other items:

  -------------------------------------------------------------------------------
  Subdirectory                        Examples of contents and use
  ----------------------------------- -------------------------------------------
  `1_Mutation_profiling/`             Position-level and region-level mutation
                                      tables, including `Pos_mut_thresh_0.03.tsv`
                                      and `Pos_mut_thresh_0.03_region.tsv`, as
                                      well as gene-wise summaries.

  `2_iclip/`                          Processed iCLIP peak data, including
                                      `SRSF3_merged.3nt_3nt.peaks.bed.gz`.

  `3_Isoform_counts/`                 Transcript-level read-count tables for the
                                      documented samples and replicates,
                                      including
                                      `filt_q8_read_counts_per_transcript.tsv`.

  `4_Co-occurence_sites/`             Read-level mutation-location and
                                      co-occurrence tables, including
                                      `filt_q8_mut_loc.tsv.gz` and
                                      `df_sites_per_read.tsv.gz`.
  -------------------------------------------------------------------------------

The exact sample, reference, replicate, and filtering labels are encoded
in the nested directory names. Retaining this directory structure is
important because the figure scripts refer to these analysis-specific
inputs.

## Software and dependencies

The workflows use a combination of command-line bioinformatics tools,
Python, and R. Requirements vary by analysis.

**Tools and packages explicitly documented in the workflow READMEs
include:**

-   Oxford Nanopore Dorado or Guppy (basecalling).
-   minimap2 (read alignment).
-   samtools, bedtools, and pigz (alignment processing).
-   R and the Bioconductor package Bambu (isoform assignment).
-   gffread (transcript sequence extraction).
-   `perbase` (per-base profiling).
-   Python packages including pandas, numpy, pysam, and biopython.
-   R packages including GenomicRanges, ggplot2, data.table, dplyr,
    matrixStats, optparse, ggupset, ComplexUpset, and patchwork, as
    required by individual analyses.

This is a folder-level summary, not a single environment
specification. Please use the README in each workflow directory to
install the packages and tools required for that particular step.
Software versions and parameters should be recorded for any independent
rerun.

## Reference genome and annotation

The workflow documentation includes examples using the mouse mm9
reference genome and Gencode annotation files (including `gencode.vM1`
resources). These reference files are not included in the supplied
archive. Users reproducing the analyses should obtain the same reference
genome and annotation release used for the relevant analysis and ensure
that chromosome naming and transcript identifiers are consistent across
FASTA, GTF, and alignment files.

Check the relevant workflow README for the exact reference-file names
and paths expected by each script.

## How to use this repository

1.  **Start with the relevant workflow README.** Each processing
    directory contains a more detailed README with command-line examples
    and expected inputs/outputs.
2.  **Use the example FASTQ files for an initial test.** The files in
    `1_Test_data/` are subsampled examples and are useful for checking
    that the relevant tools and scripts run.
3.  **Prepare the required references and input files.** Full processing
    requires the appropriate sequencing data, reference genome,
    annotation files, and software dependencies.
4.  **Run the processing workflows in sequence where needed.** Genome
    mapping precedes the genome-based downstream analyses; isoform
    mapping is a separate transcript-assignment and alignment workflow.
5.  **Generate manuscript figures.** Use the scripts in
    `figure_scripts/` and the corresponding processed inputs in
    `figure_scripts/figure_data/`. Read the README beside each figure
    script before running it.

The example commands in the individual README files may use
project-specific directory layouts (for example, `Result/` and
`genome/`). Adjust paths to match your local working directory while
keeping the expected input/output relationships unchanged.

## Reproducibility notes

-   The scripts include analysis-specific parameters and sample
    definitions. Refer to the individual README files rather than
    assuming that one set of thresholds applies to every analysis.
-   The mutation-profiling workflows distinguish analyses with and
    without deletion inclusion.
-   The isoform workflow uses uniquely assigned reads for
    isoform-specific FASTQ generation, as described in its
    documentation.
-   Figure-generation scripts are designed around the processed input
    files and directory structure documented alongside them.
-   The supplied test FASTQ files are subsampled; results obtained from
    them should not be interpreted as reproducing the full manuscript
    analyses.
-   For exact reproduction, retain the reference genome/annotation
    release, input data, software versions, command-line parameters, and
    random subsampling settings used in the original analysis.
