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
