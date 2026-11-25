#!/usr/bin/env python3
"""
Read-wise Mutation Site Co-occurrence Analysis

Purpose:
    Extract mutation positions from transcript-aligned reads and
    generate per-read mutation-location tables for co-occurrence analyses.

Preserve the original manuscript-specific implementation and parameters.
"""

############################################################
# Extract mutation locations per read from BAM alignments
# and write mutation-location tables for co-occurrence
# analysis
############################################################

##############################
# Activate environment
##############################
# conda activate isoquant

##############################
# Import required packages
##############################
import pysam
import pandas as pd
import numpy as np

from collections import defaultdict
import gzip
import subprocess
import os

############################################################
# Function: Create output directory if it does not exist
############################################################
def create_output_folder(directory):

    try:
        if not os.path.exists(directory):
            os.makedirs(directory)

    except OSError:
        print("Error: Creating directory. " + directory)


############################################################
# Function: Extract mutation positions from a single read
############################################################
def get_mutation_positions_per_read(read):

    ########################################################
    # Basic read information
    ########################################################
    reference_name = read.reference_name
    read_id = read.query_name

    ########################################################
    # Obtain aligned read/reference positions
    ########################################################
    aligned_pairs = read.get_aligned_pairs(with_seq=True)

    alignment_df = pd.DataFrame(
        aligned_pairs,
        columns=["read_pos", "ref_pos", "base"]
    )

    alignment_df = alignment_df.dropna(
        subset=["ref_pos", "base"]
    )

    ########################################################
    # Deletions:
    # read_pos = NA
    #
    # Mismatches:
    # reference base is lowercase in MD-tag output
    ########################################################
    mutation_df = alignment_df[
        alignment_df["read_pos"].isna() |
        alignment_df["base"].str.islower()
    ]

    ########################################################
    # Store mutation locations
    ########################################################
    output_df = pd.DataFrame({
        "read_id": read_id,
        "ref_name": reference_name,

        # Convert to 1-based coordinates
        "ref_mod_pos":
            mutation_df["ref_pos"].astype(int) + 1
    }).reset_index(drop=True)

    ########################################################
    # Collapse mutation positions per read
    ########################################################
    merged_df = (
        output_df
        .groupby(
            ["read_id", "ref_name"]
        )["ref_mod_pos"]
        .apply(
            lambda positions:
            ",".join(map(str, positions))
        )
        .reset_index(
            name="ref_mod_pos"
        )
    )

    return merged_df


############################################################
# Function: Extract mutation positions from all reads
############################################################
def perform_cooccurrence_analysis(
        result_directories,
        bam_filename,
        batch_size=100000,
        threads=2):
    """
    Stream reads from BAM files and write mutation
    locations to a gzipped TSV file.

    Data are written in batches to minimize memory usage.
    """

    md_bam_filename = bam_filename.replace(
        ".bam",
        "_md.bam"
    )

    for result_directory in result_directories:

        create_output_folder(result_directory)

        bam_path = os.path.join(
            result_directory,
            md_bam_filename
        )

        output_path = os.path.join(
            result_directory,
            md_bam_filename.replace(
                "_md.bam",
                "_mut_loc.tsv.gz"
            )
        )

        ####################################################
        # Open BAM file
        ####################################################
        bam_file = pysam.AlignmentFile(
            bam_path,
            "rb",
            threads=threads
        )

        ####################################################
        # Write compressed output
        ####################################################
        with gzip.open(output_path, "wt") as gz_file:

            mutation_dfs = []

            read_counter = 0
            header_written = False

            for read in bam_file.fetch(
                    until_eof=True):

                if read.is_unmapped:
                    continue

                ################################################
                # Extract mutation positions for read
                ################################################
                read_mutation_df = (
                    get_mutation_positions_per_read(
                        read
                    )
                )

                if (
                    read_mutation_df is not None and
                    not read_mutation_df.empty
                ):
                    mutation_dfs.append(
                        read_mutation_df
                    )

                read_counter += 1

                ################################################
                # Flush batch to disk
                ################################################
                if read_counter % batch_size == 0:

                    if mutation_dfs:

                        batch_df = pd.concat(
                            mutation_dfs,
                            axis=0,
                            ignore_index=True
                        )

                        batch_df.to_csv(
                            gz_file,
                            sep="\t",
                            index=False,
                            header=(
                                not header_written
                            )
                        )

                        header_written = True
                        mutation_dfs.clear()

            ####################################################
            # Write remaining data
            ####################################################
            if mutation_dfs:

                batch_df = pd.concat(
                    mutation_dfs,
                    axis=0,
                    ignore_index=True
                )

                batch_df.to_csv(
                    gz_file,
                    sep="\t",
                    index=False,
                    header=(
                        not header_written
                    )
                )

                mutation_dfs.clear()

        bam_file.close()


############################################################
# Function: Add MD and NM tags to BAM files
############################################################
def add_md_tags(
        bam_directories,
        result_directories,
        bam_filename,
        reference_fasta):

    for bam_directory, result_directory in zip(
            bam_directories,
            result_directories):

        create_output_folder(
            result_directory
        )

        input_bam = os.path.join(
            bam_directory,
            bam_filename
        )

        output_bam = os.path.join(
            result_directory,
            bam_filename.replace(
                ".bam",
                "_md.bam"
            )
        )

        ####################################################
        # Add MD and NM tags
        ####################################################
        command = (
            "samtools calmd -b "
            + input_bam
            + " "
            + reference_fasta
            + " > "
            + output_bam
        )

        print(command)

        subprocess.call(
            command,
            shell=True
        )

        ####################################################
        # Index BAM file
        ####################################################
        command = (
            "samtools index "
            + output_bam
        )

        print(command)

        subprocess.call(
            command,
            shell=True
        )


############################################################
# Project settings
############################################################
project_directory = (
    "/data/groups/anko/mohamed/"
    "oxford_nanopore/"
    "Submission/Tools/"
)

############################################################
# Input BAM files
############################################################
bam_directories = [
    "Result/4_Isoform_Mapping/"
    "SRSF3_dcDNA/rep1_1pct/"
]

bam_filename = "filt_q8.bam"

############################################################
# Output directories
############################################################
result_directories = [
    "Result/5_Co-occurence_sites/"
    "3_SRSF3_DIRCLIP/"
    "Co-SRSF3_dcDNA/rep1_1pct/"
]

############################################################
# Reference annotation
############################################################
gtf = "gencode"

reference_fasta = (
    "genome/mm9/Annotation/"
    "gencode.vM1.annotation_"
    "transcripts_with_chr.fa"
)

############################################################
# Change to project directory
############################################################
os.chdir(project_directory)

############################################################
# Step 1: Add MD/NM tags
############################################################
add_md_tags(
    bam_directories,
    result_directories,
    bam_filename,
    reference_fasta
)

############################################################
# Step 2: Extract mutation locations
############################################################
perform_cooccurrence_analysis(
    result_directories,
    bam_filename
)

############################################################
# End of script
############################################################