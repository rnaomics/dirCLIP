#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read length statistics from BAM files.

Generates:
- read_length.txt
- Read_length.xlsx

Reports:
- Mean read length
- Median read length
- Number of mapped reads
"""

import os
import subprocess
import itertools
import statistics
import pandas as pd


def create_folder(directory):
    os.makedirs(directory, exist_ok=True)


def get_read_length(read_length_file):
    """Calculate mean, median and read count from samtools length summary."""

    try:
        data = pd.read_csv(read_length_file, sep="\t", header=None)
        data.columns = ['x']

        lengths = []
        counts = []

        for entry in data['x']:
            fields = entry.split(' ')
            lengths.append(int(fields[-1]))
            counts.append(int(fields[-2]))

        expanded = [[lengths[i]] * counts[i] for i in range(len(lengths))]
        read_lengths = list(itertools.chain(*expanded))

        mean_rl = round(statistics.mean(read_lengths))
        median_rl = round(statistics.median(read_lengths))
        n_reads = sum(counts)

    except Exception:
        mean_rl = 0
        median_rl = 0
        n_reads = 0

    return mean_rl, median_rl, n_reads


def get_bam_read_length(bam_files):
    """Generate read length statistics for BAM files."""

    for bam_file in bam_files:

        directory = os.path.dirname(bam_file)

        read_length_file = os.path.join(directory, 'read_length.txt')
        excel_file = os.path.join(directory, 'Read_length.xlsx')

        command = (
            'samtools view -F 2308 ' + bam_file +
            " | cut -f10 | awk '{print length}' | sort -n | uniq -c > " +
            read_length_file
        )

        print(command)
        subprocess.call(command, shell=True)

        mean_rl, median_rl, n_reads = get_read_length(read_length_file)

        df = pd.DataFrame({
            'Mean_RL': [mean_rl],
            'Median_RL': [median_rl],
            'N_reads': [n_reads]
        })

        df.to_excel(excel_file, index=False)

        print(f'Saved: {excel_file}')


PROJECT_DIR = '/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/'

BAM_FILES = [
    'Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/sorted.bam',
    'Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/sorted.bam',
    'Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/sorted.bam',
    'Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/sorted.bam',
    'Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/sorted.bam',
    'Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/sorted.bam'
]

if __name__ == '__main__':
    os.chdir(PROJECT_DIR)
    get_bam_read_length(BAM_FILES)
