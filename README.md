<<<<<<< HEAD
# dirCLIP: Computational workflows and analysis code

This repository contains the computational workflows, analysis scripts,
and selected supporting data for the direct RNA (dRNA) and direct cDNA
(dcDNA) sequencing analyses associated with the dirCLIP study.

![Overview of the dirCLIP direct RNA and direct cDNA computational
workflows](dirCLIP_workflow_figure.png)

*Figure 1. Overview of the direct RNA (dRNA; left) and direct cDNA
(dcDNA; right) computational workflows. This is a high-level schematic:
the dcDNA basecalling directory contains both Dorado and Guppy wrappers,
and the basecaller used for each dRNA dataset is specified in the
dRNAseq documentation. See the component READMEs for exact scripts,
inputs, parameters, and software versions.*

The repository is organised into two main analysis components:

-   **[`dRNAseq/`](dRNAseq/)** -- direct RNA sequencing analysis,
    including the published analysis pipeline, the 2026 reanalysis, and
    supporting analyses and validation checks.
-   **[`dcDNAseq/`](dcDNAseq/)** -- direct cDNA sequencing analysis,
    including basecalling, genome mapping, isoform mapping, mutation
    profiling, and downstream figure-generation analyses.

This top-level README is a guide to the repository. Detailed
instructions are provided in the README files within the relevant
component and analysis directories.

> **Manuscript citation:** *dirCLIP PROFILES VARIANT-SPECIFIC
> RNA-PROTEIN INTERACTIONS VIA NANOPORE LONG-READ SEQUENCING*\
> **DOI:** To be added when available.

------------------------------------------------------------------------

## Repository structure

``` text
.
├── README.md
├── LICENSE
├── dirCLIP_workflow_figure.png
│
├── dcDNAseq/
│   ├── README.md
│   ├── 1_Test_data/
│   ├── 2_Basecalling/
│   ├── 3_Genome_Mapping/
│   ├── 4_Isoform_Mapping/
│   ├── 5_Mutation_profiling/
│   └── figure_scripts/
│       ├── 1_Read_length/
│       ├── 2_Mutation_profiling/
│       ├── 3_iCLIP_Motif_analysis/
│       ├── 4_dirCLIP_Motif_analysis/
│       ├── 5_Isoforms_detected/
│       ├── 6_Co-occurence_analysis/
│       └── figure_data/
│
└── dRNAseq/
    ├── README.md
    ├── MANIFEST.tsv
    ├── pipeline/
    ├── reanalysis_2026/
    ├── iclip_peaks/
    ├── iclip_overlap_published_sites/
    ├── supplementary_table1_checks/
    ├── nanocompore_pvalue/
    ├── pore_usage/
    └── depth_table/
```

## 1. Direct cDNA sequencing analysis (`dcDNAseq/`)

**Start here:** [`dcDNAseq/README.md`](dcDNAseq/README.md)

The `dcDNAseq/` directory contains the direct cDNA sequencing workflows
and downstream analyses.

  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Directory                                                   Contents                Main scripts / documentation
  ----------------------------------------------------------- ----------------------- ---------------------------------------------------------------------------------------
  [`1_Test_data/`](dcDNAseq/1_Test_data/)                     Subsampled example      [`README_1_Test_data.md`](dcDNAseq/1_Test_data/README_1_Test_data.md)
                                                              FASTQ files: three      
                                                              SRSF3 dcDNA replicates  
                                                              (1%) and three TotalRNA 
                                                              dcDNA no-UV replicates  
                                                              (0.2%).                 

  [`2_Basecalling/`](dcDNAseq/2_Basecalling/)                 Basecalling of Oxford   `dorado_basecall/dorado_basecaller.py`; `guppy_basecall/guppy_basecall.py`
                                                              Nanopore raw signal     
                                                              data using Dorado or    
                                                              Guppy.                  

  [`3_Genome_Mapping/`](dcDNAseq/3_Genome_Mapping/)           Genome alignment, BAM   `align_minimap.py`; `bam_merge_and_subsample.py`; see the two README files in this
                                                              processing, replicate   directory.
                                                              merging, and            
                                                              subsampling.            

  [`4_Isoform_Mapping/`](dcDNAseq/4_Isoform_Mapping/)         Bambu-based             `isoform_pipeline.py`; `bambu_for_read_assignment.R`;
                                                              read-to-transcript      [`README_isoform_pipeline.md`](dcDNAseq/4_Isoform_Mapping/README_isoform_pipeline.md)
                                                              assignment,             
                                                              isoform-specific read   
                                                              processing/alignment,   
                                                              and transcript-level    
                                                              read counting.          

  [`5_Mutation_profiling/`](dcDNAseq/5_Mutation_profiling/)   Per-base mutation       `mutation_profiling_pipeline.py`; `significant_mutations.R`; see the two README files
                                                              profiling and           in this directory.
                                                              statistical             
                                                              identification of       
                                                              significant mutation    
                                                              sites.                  

  [`figure_scripts/`](dcDNAseq/figure_scripts/)               R and Python scripts    See the analysis-specific README files in each subdirectory.
                                                              for downstream analyses 
                                                              and manuscript figures. 
  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Figure-generation and downstream analysis directories

  ----------------------------------------------------------------------------------------------------------------------
  Directory                                                                          Analysis
  ---------------------------------------------------------------------------------- -----------------------------------
  [`1_Read_length/`](dcDNAseq/figure_scripts/1_Read_length/)                         Read-length statistics and
                                                                                     read-length distribution plots.

  [`2_Mutation_profiling/`](dcDNAseq/figure_scripts/2_Mutation_profiling/)           Mutation-rate comparisons,
                                                                                     mutation-site overlap, gene-level
                                                                                     summaries, iCLIP overlap,
                                                                                     exon-junction metaplots,
                                                                                     transcript-position distributions,
                                                                                     mutation-site annotation, and
                                                                                     dcDNA/dRNA overlap.

  [`3_iCLIP_Motif_analysis/`](dcDNAseq/figure_scripts/3_iCLIP_Motif_analysis/)       Sequence extraction around iCLIP
                                                                                     sites and motif-enrichment
                                                                                     analysis.

  [`4_dirCLIP_Motif_analysis/`](dcDNAseq/figure_scripts/4_dirCLIP_Motif_analysis/)   Crosslink-site and
                                                                                     background-sequence extraction,
                                                                                     motif enrichment, and CNNC
                                                                                     motif-distance analysis.

  [`5_Isoforms_detected/`](dcDNAseq/figure_scripts/5_Isoforms_detected/)             Quantification of detected
                                                                                     transcript isoforms per gene.

  [`6_Co-occurence_analysis/`](dcDNAseq/figure_scripts/6_Co-occurence_analysis/)     Read-level mutation-site
                                                                                     co-occurrence and inter-site
                                                                                     distance analyses.

  [`figure_data/`](dcDNAseq/figure_scripts/figure_data/)                             Selected processed tables used by
                                                                                     downstream figure scripts.
  ----------------------------------------------------------------------------------------------------------------------

The figure-analysis directories contain README files describing the
corresponding scripts and their inputs/outputs. Not every script has a
separate one-to-one README; some workflows are documented by a README
covering multiple scripts.

## 2. Direct RNA sequencing analysis (`dRNAseq/`)

**Start here:** [`dRNAseq/README.md`](dRNAseq/README.md)

The `dRNAseq/` directory contains the published direct RNA analysis
workflow, the 2026 reanalysis, and supporting analyses and validation
checks.

  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Directory                                                                    Contents                     Documentation
  ---------------------------------------------------------------------------- ---------------------------- ----------------------------------------------------------------------------------------------
  [`pipeline/`](dRNAseq/pipeline/)                                             Published analysis pipeline  [`pipeline/README.md`](dRNAseq/pipeline/README.md)
                                                                               (scripts 1--10), covering    
                                                                               basecalling, mapping,        
                                                                               Nanopolish/NanopolishComp,   
                                                                               Nanocompore signal-deviation 
                                                                               analysis, read-termination   
                                                                               analysis, and merging of the 
                                                                               two channels.                

  [`reanalysis_2026/`](dRNAseq/reanalysis_2026/)                               Slurm scripts for the 2026   [`reanalysis_2026/README.md`](dRNAseq/reanalysis_2026/README.md)
                                                                               signal-level reanalysis of   
                                                                               all nine libraries,          
                                                                               including mapping,           
                                                                               eventalign, NanopolishComp,  
                                                                               Nanocompore, and             
                                                                               significant-site extraction. 

  [`iclip_peaks/`](dRNAseq/iclip_peaks/)                                       Processing of published      [`iclip_peaks/README.md`](dRNAseq/iclip_peaks/README.md)
                                                                               SRSF3 iCLIP peaks, liftOver  
                                                                               from mm9 to GRCm39,          
                                                                               peak-confidence tiers, and   
                                                                               read-depth lookup            
                                                                               generation.                  

  [`iclip_overlap_published_sites/`](dRNAseq/iclip_overlap_published_sites/)   Overlap analyses between     [`iclip_overlap_published_sites/README.md`](dRNAseq/iclip_overlap_published_sites/README.md)
                                                                               published site sets and      
                                                                               SRSF3 iCLIP peaks, including 
                                                                               odds ratios, bootstrap       
                                                                               analyses, recovery           
                                                                               statistics, and              
                                                                               supplementary-table          
                                                                               generation.                  

  [`supplementary_table1_checks/`](dRNAseq/supplementary_table1_checks/)       Extraction of Supplementary  [`supplementary_table1_checks/README.md`](dRNAseq/supplementary_table1_checks/README.md)
                                                                               Table 1 and checks of region 
                                                                               definitions and              
                                                                               cross-channel merging rules. 

  [`nanocompore_pvalue/`](dRNAseq/nanocompore_pvalue/)                         Checks and source evidence   [`nanocompore_pvalue/README.md`](dRNAseq/nanocompore_pvalue/README.md)
                                                                               concerning                   
                                                                               Benjamini--Hochberg          
                                                                               adjustment of Nanocompore    
                                                                               1.0.x `GMM_logit_pvalue`     
                                                                               values.                      

  [`pore_usage/`](dRNAseq/pore_usage/)                                         Processing of MinKNOW run    [`pore_usage/README.md`](dRNAseq/pore_usage/README.md)
                                                                               reports and generation of    
                                                                               pore-usage source data and   
                                                                               supplementary plots.         

  [`depth_table/`](dRNAseq/depth_table/)                                       FASTQ statistics and         [`depth_table/README.md`](dRNAseq/depth_table/README.md)
                                                                               generation of the            
                                                                               supplementary                
                                                                               sequencing-depth and         
                                                                               site-yield table.            

  [`MANIFEST.tsv`](dRNAseq/MANIFEST.tsv)                                       File inventory with file     ---
                                                                               sizes and SHA-256 checksums  
                                                                               for the contents of          
                                                                               `dRNAseq/` (excluding the    
                                                                               manifest itself).            
  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Important distinction between the dRNA workflows

The published workflow in `dRNAseq/pipeline/` and the 2026 reanalysis in
`dRNAseq/reanalysis_2026/` are related but distinct analysis paths. The
dRNAseq README describes which scripts and basecalled data belong to
each path. Please consult that README before attempting to reproduce the
reported results.

------------------------------------------------------------------------

## 3. Where to find a particular analysis

  ------------------------------------------------------------------------------
  If you want to find...              Start here
  ----------------------------------- ------------------------------------------
  dcDNA basecalling                   `dcDNAseq/2_Basecalling/`

  dcDNA genome alignment and BAM      `dcDNAseq/3_Genome_Mapping/`
  processing                          

  dcDNA isoform assignment and        `dcDNAseq/4_Isoform_Mapping/`
  counting                            

  dcDNA mutation-site calling         `dcDNAseq/5_Mutation_profiling/`

  dcDNA manuscript figure scripts     `dcDNAseq/figure_scripts/` and the README
                                      in the relevant analysis directory

  Processed inputs for dcDNA figure   `dcDNAseq/figure_scripts/figure_data/`
  analyses                            

  Published direct RNA pipeline       `dRNAseq/pipeline/`

  2026 direct RNA reanalysis          `dRNAseq/reanalysis_2026/`

  iCLIP peak processing and overlap   `dRNAseq/iclip_peaks/` and
  analyses                            `dRNAseq/iclip_overlap_published_sites/`

  Supplementary Table 1 checks        `dRNAseq/supplementary_table1_checks/`

  Nanocompore p-value adjustment      `dRNAseq/nanocompore_pvalue/`
  checks                              

  Pore-usage analysis                 `dRNAseq/pore_usage/`

  Sequencing-depth table              `dRNAseq/depth_table/`
  ------------------------------------------------------------------------------

------------------------------------------------------------------------

## 4. Data and reference resources

This repository contains analysis scripts, selected example data, and
selected processed tables. It is **not** a complete archive of all raw
sequencing data or external reference resources.

### Included

-   Subsampled example FASTQ files in `dcDNAseq/1_Test_data/`.
-   Selected processed tables used by dcDNA figure-generation scripts in
    `dcDNAseq/figure_scripts/figure_data/`.
-   Supporting source data and analysis outputs explicitly present in
    the dRNAseq directories.
-   README files documenting workflow-specific inputs, commands, and
    outputs.

### Not included

-   Full raw sequencing runs (FAST5/POD5 and complete FASTQ datasets).
-   All reference genomes, transcriptomes, annotations, and chain files
    required by the workflows.
-   The external MinKNOW run reports and other project-specific input
    files described in the component READMEs.
-   The processed dRNA site tables identified in the dRNAseq
    documentation as manuscript supplementary tables.

The dcDNA `figure_data/` directory contains selected processed inputs,
but it does not replace every external input needed by every figure
script. For example, some analyses require BAM files, transcript
annotations, or other intermediate files generated by upstream
workflows. Check the README for the specific analysis before assuming
that a figure can be regenerated from the files in `figure_data/` alone.

**Raw data availability:** Add the relevant public repository/accession
numbers here when available.

------------------------------------------------------------------------

## 5. Software and computational environment

The scripts use a combination of Python, R, shell, and Slurm.
Dependencies vary by workflow; there is no single environment file that
installs all requirements.

Examples of tools used across the repository include:

-   **Basecalling:** Oxford Nanopore Dorado and Guppy.
-   **Alignment and BAM processing:** minimap2, samtools, bedtools, and
    pigz.
-   **dcDNA analysis:** Bambu, gffread, `perbase`, Python packages
    (including pandas, numpy, pysam, and biopython), and R/Bioconductor
    packages (including GenomicRanges, ggplot2, data.table, and dplyr).
-   **dRNA analysis:** Nanopolish, NanopolishComp, Nanocompore, R2Dtool,
    UCSC liftOver, Python packages (including numpy, pandas, scipy,
    statsmodels, openpyxl, and matplotlib), and R packages (including
    tidyverse and glmGamPoi).

For exact versions, installation instructions, and workflow-specific
dependencies, consult the relevant component README:

-   [`dcDNAseq/README.md`](dcDNAseq/README.md)
-   [`dRNAseq/README.md`](dRNAseq/README.md)

The dRNAseq README records the versions used for the 2026 reanalysis.
The scripts in the published dRNA pipeline do not all pin software
versions.

------------------------------------------------------------------------

## 6. Reproducing an analysis

1.  Identify whether the analysis belongs to the direct cDNA
    (`dcDNAseq/`) or direct RNA (`dRNAseq/`) component.
2.  Read the component README.
3.  Open the README in the relevant workflow or figure-analysis
    directory.
4.  Obtain the required sequencing data, reference files, and other
    external inputs listed in that README.
5.  Install the software versions and dependencies required for that
    workflow.
6.  Run the scripts in the documented order, using paths appropriate for
    your computing environment.

The example FASTQ files in `dcDNAseq/1_Test_data/` are subsampled and
are intended for testing selected steps. They should not be expected to
reproduce the full manuscript results.

Some dcDNA scripts retain paths or working-directory assumptions from
the original analysis environment. These may need to be adapted before
running the scripts on another system; check each script's configuration
section and its local workflow README.

------------------------------------------------------------------------

## 7. File and folder naming

Some historical directory and script names contain spelling variants
(for example, `6_Co-occurence_analysis/` and `4_Co-occurence_sites/`).
These names are retained to preserve compatibility with existing scripts
and documentation.

------------------------------------------------------------------------

## Citation

Please cite the associated dirCLIP manuscript when using the methods or
code in this repository:

**dirCLIP PROFILES VARIANT-SPECIFIC RNA-PROTEIN INTERACTIONS VIA
NANOPORE LONG-READ SEQUENCING**

**DOI:** To be added when available.

## License

This project is distributed under the MIT License. See
[`LICENSE`](LICENSE).

## Contact

For questions about the code or analyses, please contact the
corresponding author or repository maintainer listed with the
manuscript.
=======
# dirCLIP

**dirCLIP** is a computational toolkit designed for analyzing and profiling RNA-protein interactions using nanopore sequencing data. This repository contains scripts and workflows used in the manuscript submitted to *Nature Methods*.

## 📁 Folder Overview

- `Alignment/`
  - Contains scripts and tools for aligning nanopore reads to reference genomes.

- `Motif_analysis/`
  - Includes R scripts for analyzing sequence motifs and CNNC distances from mutation sites.

- `Mutation_profiling/`
  - Provides Python and R scripts for identifying, grouping, and visualizing mutation sites across samples.

## 🚀 Usage Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/rnaomics/dirCLIP.git
   cd dirCLIP
   ```

2. Navigate to the desired folder and run the scripts using Python or R:

To find dirCLIP sites using Mutation profiling, go to the folder Mutation profiling.
   ```bash
   cd Mutation_profiling
   ```

3. Ensure required dependencies are installed (e.g., Bioconductor packages for R, pandas/numpy for Python).

## 📜 License
This project is released under the MIT License.

## 📬 Contact
For questions, please contact Mohamed Bahrudeen (mailto:mohamed.mohamedbahrdueen@tuni.fi).
>>>>>>> 90ef3501dc4c2b0f9aeb9739ba7950194de15365
