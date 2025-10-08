# -*- coding: utf-8 -*-
"""
bam_processing_pipeline.py

Pipeline for merging BAM files, subsampling reads, and converting to BED format.
Designed for Oxford Nanopore sequencing data using samtools and bedtools.

Author: Mohamed B
Created: June 15, 2023
"""

import os
import subprocess

def create_folder(directory):
    """
    Creates a directory if it does not exist.

    Args:
        directory (str): Path to the target directory.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def merge_bam_files(Bam_dirs, Result_dirs, bam_file, threads):
    """
    Merges BAM files from multiple replicate directories into a single sorted BAM.

    Args:
        Bam_dirs (list): List of lists containing replicate BAM directory paths.
        Result_dirs (list): List of output directories for merged BAM files.
        bam_file (str): Name of the BAM file to merge.
        threads (int): Number of threads for samtools operations.
    """
    for i, Bam_dir_group in enumerate(Bam_dirs):
        bam_files = [os.path.join(dir_path, bam_file) for dir_path in Bam_dir_group]
        result_dir = Result_dirs[i]
        merge_file = os.path.join(result_dir, bam_file)
        create_folder(result_dir)

        # Merge BAM files
        command_merge = f"samtools merge -f -@ {threads} {merge_file} {' '.join(bam_files)}"
        subprocess.call(command_merge, shell=True)

        # Remove temporary sort files
        command_rm = f"rm {merge_file}.tmp.*"
        subprocess.call(command_rm, shell=True)

        # Sort the merged BAM file
        command_sort = f"samtools sort -@ {threads} -o {merge_file} {merge_file}"
        subprocess.call(command_sort, shell=True)

def subsample_bam(N, sample_percent, Result_dirs, bam_file, bed_file, threads):
    """
    Subsamples BAM files and converts them to BED format.

    Args:
        N (int): Number of subsamples to generate.
        sample_percent (int): Percent of reads to retain in each subsample.
        Result_dirs (list): Directories containing merged BAM files.
        bam_file (str): Name of the BAM file to subsample.
        bed_file (str): Output BED file name.
        threads (int): Number of threads for sorting.
    """
    seeds = list(range(100, N * 100 + 1, 100))  # Unique seeds for each subsample

    for result_dir in Result_dirs:
        for i in range(N):
            subsample_dir = os.path.join(result_dir, f"subsample{i+1}")
            create_folder(subsample_dir)

            merge_file = os.path.join(result_dir, bam_file)
            subsample_file = os.path.join(subsample_dir, bam_file)
            subsample_bed_file = os.path.join(subsample_dir, bed_file)

            # Subsample and sort BAM
            command_view_sort = f"samtools view -s {seeds[i]}.{sample_percent} -b {merge_file} | samtools sort -@ {threads} -o {subsample_file}"
            subprocess.call(command_view_sort, shell=True)

            # Index BAM
            command_index = f"samtools index {subsample_file}"
            subprocess.call(command_index, shell=True)

            # Convert BAM to BED and compress
            command_bed = f"bedtools bamtobed -i {subsample_file} | pigz > {subsample_bed_file}"
            subprocess.call(command_bed, shell=True)

# =================== USER CONFIGURATION ===================

# Set working project directory
Project_dir = "/data/groups/anko/mohamed/oxford_nanopore/"
os.chdir(Project_dir)

# File names
bam_file = "sorted.bam"
bed_file = "sorted.bed.gz"
threads = 2  # Number of samtools threads
N = 3        # Number of subsamples
sample_percent = 70  # % of reads to retain per subsample

# Input BAM replicate directories grouped by experiment
Bam_dirs = [
    ["Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1A_deltaC_CLIP_cDNA_rep1/",
     "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1A_deltaC_CLIP_cDNA_rep2/"],
    ["Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_rep1/",
     "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_rep2/"],
    ["Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1A_deltaD_CLIP_cDNA_rep1/"]
]

# Output directories for merged BAM files
Result_dirs = [
    "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1A_deltaC_CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/Hmga1A_deltaD_CLIP_cDNA_merge/"
]

# Run pipeline
merge_bam_files(Bam_dirs, Result_dirs, bam_file, threads)
subsample_bam(N, sample_percent, Result_dirs, bam_file, bed_file, threads)
