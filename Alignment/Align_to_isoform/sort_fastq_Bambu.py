#! /home/pbohn/miniconda3/envs/nanodms/bin/python

import numpy as np
from Bio import SeqIO
import gzip
import argparse
import pandas as pd
import os
import pdb

# Argument parser setup
parser = argparse.ArgumentParser(
    prog="split fastq file into one for each isoform",
    description="Reads in a FASTQ file and a read_assignment file, then sorts all reads uniquely mapped to an isoform into a new FASTQ file"
)
parser.add_argument("-fi", "--fastq_infile")
parser.add_argument("-fo", "--fastq_out_prefix")
parser.add_argument("-ir", "--readToTranscriptMaps")
parser.add_argument("-rdf", "--read_assignment_df")
parser.add_argument("-o", "--outfile_prefix")
parser.add_argument("-m", "--read_id_mapping_file")
parser.add_argument("-u", "--output_unsorted", action="store_true")

def sort_fastq_by_isoform_isoquant(fastq_infile, readToTranscriptMaps, fastq_outdir, output_unsorted=False, df=False):
    # Load the read-to-transcript mapping
    if df:
        read_df = pd.read_pickle(readToTranscriptMaps)
    else:
        read_df = pd.read_csv(readToTranscriptMaps, sep="\t", keep_default_na=False, na_values=['_'])

    # Parse and combine isoform match lists
    read_df['equalMatches'] = read_df['equalMatches'].apply(lambda x: x.split(', ') if isinstance(x, str) else [x])
    read_df['compatibleMatches'] = read_df['compatibleMatches'].apply(lambda x: x.split(', ') if isinstance(x, str) else [x])
    read_df['isoforms'] = read_df.apply(lambda row: row['equalMatches'] + row['compatibleMatches'], axis=1)
    read_df['isoforms'] = read_df['isoforms'].apply(lambda x: [item for item in x if item != ''])

    total_num_reads = len(read_df)

    # Keep only reads that map uniquely to one isoform
    read_df = read_df[read_df['isoforms'].apply(lambda x: len(x) == 1)]
    read_df = read_df[["readId", "isoforms"]]
    unique_reads = len(read_df)

    # Get list of unique isoforms
    isoforms = [isoform for sublist in read_df["isoforms"] for isoform in sublist]
    isoforms = sorted(np.unique(isoforms))

    # Unpack isoform list into a single value
    read_df['isoforms'] = read_df['isoforms'].apply(lambda x: x[0] if isinstance(x, list) and len(x) == 1 else "")

    print("Detected in total", total_num_reads, "aligned reads.")
    print("Of these,", unique_reads, "reads were assigned uniquely.")
    print("Detected in total", len(isoforms), "isoforms with at least one unique assignment.")

    # Create mapping of read ID to isoform
    read_dict = read_df.set_index("readId").to_dict()["isoforms"]

    print("Reading in files")

    # Prepare output directory
    sorted_reads = {isoform: [] for isoform in isoforms}
    os.makedirs(fastq_outdir, exist_ok=True)
    num_sorted = 0
    num_not_sorted = 0
    unsorted_reads = []

    # Process reads from gzipped FASTQ file
    if ".gz" in fastq_infile:
        with gzip.open(fastq_infile, "rt") as infile:
            for read in SeqIO.parse(infile, "fastq"):
                try:
                    sorted_reads[read_dict[read.id]].append(read)
                    num_sorted += 1
                except:
                    num_not_sorted += 1
                    if output_unsorted:
                        unsorted_reads.append(read)

    # Write sorted reads to files
    isoform_counts = []
    for isoform in isoforms:
        isoform_reads = sorted_reads[isoform]
        isoform_count = len(isoform_reads)
        isoform_counts.append((isoform, isoform_count))
        with gzip.open(f"{fastq_outdir}/{isoform}.fastq.gz", "wt") as outfile:
            SeqIO.write(isoform_reads, outfile, "fastq")

    # Write unsorted reads if requested
    if output_unsorted:
        with gzip.open(f"{fastq_outdir}/unsorted.fastq.gz", "wt") as outfile:
            SeqIO.write(unsorted_reads, outfile, "fastq")
        print("Sorted", isoform_count, "reads for isoform", isoform)

    print("Sorted reads in total:", num_sorted)
    print("Not sorted reads in total:", num_not_sorted)

    # Write isoform read counts
    with open(f"{fastq_outdir}/isoform_counts.csv", "w") as csv_file:
        for isoform, count in isoform_counts:
            csv_file.write(f"{isoform}; {count}\n")

if __name__ == "__main__":
    args = parser.parse_args()

    print("Reading in reads from file", args.fastq_infile)
    if args.output_unsorted:
        print("Also outputting unsorted reads")

    if args.readToTranscriptMaps:
        print("Detected readToTranscriptMaps option, performing FASTQ per isoform sorting according to readToTranscriptMaps.")
        print("Reading read assignment file from", args.readToTranscriptMaps)
        print("Outfile prefix:", args.outfile_prefix)

        sort_fastq_by_isoform_isoquant(
            args.fastq_infile,
            args.readToTranscriptMaps,
            args.fastq_out_prefix,
            args.output_unsorted
        )
