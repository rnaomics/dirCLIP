#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Isoform-specific Read Assignment and Alignment Pipeline

Refactored command-line workflow with subcommands:

    prepare
    assign
    sort-fastq
    align
    count
    all

NOTE:
This file contains the fully refactored workflow architecture and the
preprocessing module implementation. Existing domain-specific functions
(e.g. sort_fastq_by_isoform, align_reads_to_isoform) should be placed
in the marked section below or imported from separate modules.
"""

import argparse
import gzip
import os
import subprocess
from pathlib import Path
from Bio import SeqIO
import shutil
import pandas as pd
import numpy as np
from collections import defaultdict
import pysam
import pdb

# ============================================================
# Utility functions
# ============================================================

def create_folder(directory):
    os.makedirs(directory, exist_ok=True)

def validate_inputs(bam_files, output_dirs, genome_fasta, gtf_file, r_script):
    """Validate required files and directories."""

    if len(bam_files) != len(output_dirs):
        raise ValueError(
            'Number of BAM files must match number of output directories.'
        )

    for bam in bam_files:
        if not Path(bam).exists():
            raise FileNotFoundError(f'BAM file not found: {bam}')

    if not Path(genome_fasta).exists():
        raise FileNotFoundError(genome_fasta)

    if not Path(gtf_file).exists():
        raise FileNotFoundError(gtf_file)

    if not Path(r_script).exists():
        raise FileNotFoundError(r_script)

    for outdir in output_dirs:
        os.makedirs(outdir, exist_ok=True)

# ============================================================
# PREPARE MODULE
# ============================================================

def parse_attributes(attribute_string):
    """Parse GTF attribute column into a dictionary."""

    attributes = {}

    for attribute in attribute_string.split(';'):

        attribute = attribute.strip()

        if not attribute:
            continue

        parts = attribute.split(' ', 1)

        if len(parts) != 2:
            continue

        key = parts[0]
        value = parts[1].replace('"', '')

        attributes[key] = value

    return attributes



def add_chr_to_gtf(input_gtf, output_gtf):
    """Add Chr_ prefix to transcript IDs."""

    with gzip.open(input_gtf, 'rt') as infile, gzip.open(output_gtf, 'wt') as outfile:

        for line in infile:

            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")

            attributes = parse_attributes(fields[8])

            transcript_id = attributes.get('transcript_id')
            transcript_id = f'Chr_{transcript_id.split(".")[0]}'
            if transcript_id:
                attributes['transcript_id'] = transcript_id

                fields[8] = "; ".join([f"{key} \"{value}\"" for key, value in attributes.items()])

            outfile.write("\t".join(fields) + "\n")

def get_transcriptome_from_gtf(gtf_file, genome_fasta, output_fasta):
    """Generate transcriptome FASTA using gffread."""

    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    cmd = (
        f'gunzip -c {gtf_file} | '
        f'gffread -w {output_fasta} -g {genome_fasta} -'
    )

    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

def get_transcriptwise_fasta(input_fasta, output_dir):
    """Split transcriptome FASTA into one FASTA per transcript."""

    create_folder(output_dir)

    with open(input_fasta, "r") as handle:

        for record in SeqIO.parse(handle, "fasta"):

            output_file = os.path.join(output_dir, f'{record.id}.fa')

            with open(output_file, 'w') as out:
                SeqIO.write(record, out, "fasta")

def prepare_transcriptome(gtf, genome_fasta, transcript_dir):
    """Complete transcriptome preparation workflow."""

    gtf_with_chr = gtf.replace('.gtf.gz', '_with_chr.gtf.gz')

    transcript_fasta = gtf.replace('.gtf.gz', '_transcripts_with_chr.fa')

    add_chr_to_gtf(gtf, gtf_with_chr)
    get_transcriptome_from_gtf(gtf_with_chr, genome_fasta, transcript_fasta)
    get_transcriptwise_fasta(transcript_fasta, transcript_dir)

# ------------------------------------------------------------
# Main Bambu wrapper
# ------------------------------------------------------------

def bambu_read_assignment(
        bam_files,
        output_dirs,
        genome_fasta,
        gtf_file,
        r_script='bambu_for_read_assignment.R',
        parameter_file='parameters.txt',
        remove_rdata=True):
    """Run Bambu read assignment."""

    validate_inputs(
        bam_files,
        output_dirs,
        genome_fasta,
        gtf_file,
        r_script
    )
    print('Running Bambu assignment...')
    print('Writing parameter file...')

    with open(parameter_file, 'w') as handle:
        handle.write(f'genome_file={genome_fasta}\n')
        handle.write(f'gtf_file={gtf_file}\n')
        handle.write(f'Bam_files={",".join(bam_files)}\n')
        handle.write(f'bambu_dirs={",".join(output_dirs)}\n')

    cmd = ['R', 'CMD', 'BATCH', r_script]

    print('Running Bambu read assignment...')
    print(' '.join(cmd))

    subprocess.run(cmd, check=True)

    if remove_rdata and Path('.RData').exists():
        os.remove('.RData')

    print('Finished Bambu read assignment.')

# ------------------------------------------------------------------
# FASTQ sorting by isoform
# ------------------------------------------------------------------

def sort_fastq_by_isoform(
        fastq_file,
        read_assignment_file,
        output_dir,
        output_unsorted=False):

    os.makedirs(output_dir, exist_ok=True)

    read_df = pd.read_csv(read_assignment_file, sep = "\t", keep_default_na=False, na_values=['_'])

    #### Splitting the string into list for equalMatches and compatibleMatches
    read_df['equalMatches'] = read_df['equalMatches'].apply(lambda x: x.split(', ') if isinstance(x, str) else [])
    read_df['compatibleMatches'] = read_df['compatibleMatches'].apply(lambda x: x.split(', ') if isinstance(x, str) else [])

    #### Merging the equalMatches and compatibleMatches
    read_df['isoforms'] = read_df.apply(lambda row: row['equalMatches'] + row['compatibleMatches'], axis=1)

    # Remove empty strings ('')  and "NA" from the merged lists
    read_df['isoforms'] = read_df['isoforms'].apply(lambda x: [i for i in x if i not in ['', 'NA']])

    total_reads = len(read_df)

    # Filter rows where 'isoforms' has exactly one entry
    read_df = read_df[read_df['isoforms'].apply(lambda x: len(x) == 1)]
    read_df = read_df[["readId","isoforms"]]
    unique_reads = len(read_df)

    # Flatten the list of isoforms and convert it to a 1D array
    isoforms = [isoform for sublist in read_df["isoforms"] for isoform in sublist]
    isoforms = sorted(np.unique(isoforms)) # Use np.unique to get the unique isoforms

    print(f'Total reads: {total_reads}')
    print(f'Uniquely assigned reads: {unique_reads}')
    print(f'Isoforms detected: {len(isoforms)}')

    read_df['isoforms'] = read_df['isoforms'].apply(lambda x: x[0])
    read_dict = read_df[["readId", "isoforms"]].set_index('readId').to_dict()['isoforms']

    sorted_reads = {isoform: [] for isoform in isoforms}

    unsorted_reads = []

    with gzip.open(fastq_file, 'rt') as infile:

        for read in SeqIO.parse(infile, 'fastq'):

            if read.id in read_dict:
                sorted_reads[read_dict[read.id]].append(read)
            elif output_unsorted:
                unsorted_reads.append(read)

    isoform_counts = []

    for isoform in isoforms:

        reads = sorted_reads[isoform]

        outfile = os.path.join(output_dir, f'{isoform}.fastq.gz')

        with gzip.open(outfile, 'wt') as handle:
            SeqIO.write(reads, handle, 'fastq')

        isoform_counts.append((isoform, len(reads)))

    count_file = os.path.join(output_dir, 'isoform_counts.csv')

    with open(count_file, 'w') as handle:
        for isoform, count in isoform_counts:
            handle.write(f'{isoform}; {count}\n')

    if output_unsorted:

        unsorted_file = os.path.join(output_dir, 'unsorted.fastq.gz')

        with gzip.open(unsorted_file, 'wt') as handle:
            SeqIO.write(unsorted_reads, handle, 'fastq')


def sort_fastq_reads_per_isoform(
        fastq_files,
        output_dirs,
        assignment_file="Out/readToTranscriptMaps.tsv"
):
    print('Sorting FASTQ reads by isoform...')

    if len(fastq_files) != len(output_dirs):
        raise ValueError("FASTQ and output directory counts differ.")

    for fastq_file, output_dir in zip(fastq_files, output_dirs):

        read_assignment = os.path.join(output_dir, assignment_file)

        fastq_outdir = os.path.join(output_dir, "fastq_by_isoform")

        create_folder(fastq_outdir)

        sort_fastq_by_isoform(fastq_file, read_assignment, fastq_outdir)

        src = os.path.join(fastq_outdir, "isoform_counts.csv")

        dst = os.path.join(output_dir, "isoform_counts.csv")

        shutil.move(src, dst)

        print(f"Finished sorting: {fastq_file}")
# ------------------------------------------------------------------
# Align reads to isoform
# ------------------------------------------------------------------
def align_reads_to_isoform(
        output_dirs,
        bam_name,
        transcript_dir,
        threads=4,
        sam_flag=2308,
        delete_intermediate=False
):

    if not os.path.exists(transcript_dir):
        raise FileNotFoundError(transcript_dir)

    for output_dir in output_dirs:

        fastq_path = os.path.join(output_dir, "fastq_by_isoform")

        bam_dir = os.path.join(output_dir, "bam_by_isoform")

        create_folder(bam_dir)

        #### Load isoform count ##############
        isoform_counts_file = os.path.join(output_dir, "isoform_counts.csv")

        isoform_counts = pd.read_csv(isoform_counts_file, header=None)

        isoform_counts[["isoform", "counts"]] = isoform_counts[0].str.split(";", expand=True) # Split the data into two columns
        isoform_counts = isoform_counts.drop(columns=[0]) # Drop the original combined column
        isoform_counts['counts'] = pd.to_numeric(isoform_counts['counts'])

        ##### Perform alignment ################
        for isoform in isoform_counts["isoform"]:

            isoform = isoform.strip()

            fastq_file = os.path.join(fastq_path, f"{isoform}.fastq.gz")

            fasta_file = os.path.join(transcript_dir, f"{isoform}.fa")

            if not os.path.exists(fasta_file):

                print(
                    f"WARNING: missing "
                    f"{fasta_file}"
                )

                continue

            bam_file = os.path.join(bam_dir, f"{isoform}.bam")

            cmd = (
                f"minimap2 "
                f"-ax splice "
                f"-uf "
                f"-t {threads} "
                f"-k14 "
                f"{fasta_file} "
                f"{fastq_file} "
                f"| samtools view "
                f"-bh "
                f"-F {sam_flag} "
                f"| samtools sort "
                f"-O BAM "
                f"-o {bam_file}"
            )

            subprocess.run(cmd, shell=True, check=True)

        merged_bam = os.path.join(output_dir, "merged.bam")

        merge_bam_files_in_batches(bam_dir, merged_bam, 1000)

        final_bam = os.path.join(output_dir, bam_name)

        subprocess.run(
            f"samtools sort "
            f"-o {final_bam} "
            f"{merged_bam}",
            shell=True,
            check=True
        )

        subprocess.run(
            f"samtools index "
            f"{final_bam}",
            shell=True,
            check=True
        )

        if delete_intermediate:

            shutil.rmtree(fastq_path, ignore_errors=True)
            shutil.rmtree(bam_dir, ignore_errors=True)
            os.remove(merged_bam)

# ------------------------------------------------------------------
# BAM merging
# ------------------------------------------------------------------

def merge_bam_files_in_batches(
        bam_directory,
        output_bam,
        batch_size=500,
        clean_header=True):

    bam_files = sorted([
        os.path.join(bam_directory, f)
        for f in os.listdir(bam_directory)
        if f.endswith('.bam')
    ])

    if len(bam_files) == 0:
        raise ValueError('No BAM files found.')

    if len(bam_files) <= batch_size:

        command = ('samtools merge -f ' + output_bam + ' ' + ' '.join(bam_files))

        print(command)
        subprocess.call(command, shell=True)

        if clean_header:
            remove_bam_header(output_bam)

        return

    partial_bams = []

    for i in range(0, len(bam_files), batch_size):

        batch = bam_files[i:i + batch_size]

        partial_bam = output_bam.replace('.bam', f'_partial_{i}.bam')

        partial_bams.append(partial_bam)

        command = ('samtools merge -f ' + partial_bam + ' ' + ' '.join(batch))

        print(command)
        subprocess.call(command, shell=True)

        if clean_header:
            remove_bam_header(partial_bam)

    command = ('samtools merge -f ' + output_bam + ' ' + ' '.join(partial_bams))

    print(command)
    subprocess.call(command, shell=True)

    if clean_header:
        remove_bam_header(output_bam)

    for bam in partial_bams:
        if os.path.exists(bam):
            os.remove(bam)

# ------------------------------------------------------------------
# BAM header cleanup
# ------------------------------------------------------------------

def remove_bam_header(input_bam):

    header_file = input_bam.replace('.bam', '_header.txt')

    temp_bam = input_bam.replace('.bam', '_temp.bam')

    command = (
        f"samtools view -H {input_bam} "
        f"| grep -E '^@HD|^@SQ' "
        f"> {header_file}"
    )

    print(command)
    subprocess.call(command, shell=True)

    command = (
        f"samtools reheader "
        f"{header_file} "
        f"{input_bam} "
        f"> {temp_bam}"
    )

    print(command)
    subprocess.call(command, shell=True)

    os.replace(temp_bam, input_bam)

    if os.path.exists(header_file):
        os.remove(header_file)

# ------------------------------------------------------------------
# Read counting
# ------------------------------------------------------------------

def count_reads_per_transcript(
        output_dirs,
        bam_name,
        countby,
        gtf_file):

    # --- Step 1: Parse GTF to map transcripts to gene names ---
    transcript_to_gene = {}

    with gzip.open(gtf_file, 'rt') as gtf:

        for line in gtf:

            if line.startswith('#'):
                continue

            fields = line.rstrip().split('\t')

            if fields[2] != 'transcript':
                continue

            attributes = fields[8]
            attrs = {kv.strip().split(' ')[0]: kv.strip().split(' ')[1].strip('"') 
                    for kv in attributes.strip(';').split('; ') if ' ' in kv}

            transcript_id = attrs.get('transcript_id')
            
            gene = attrs.get(countby, 'NA')

            if transcript_id:
                transcript_to_gene[transcript_id] = gene

    # --- Step 2: Count reads per transcript ---
    for output_dir in output_dirs:

        bam_file = os.path.join(output_dir, bam_name)

        output_file = os.path.join(output_dir, bam_name.replace('.bam', '_read_counts_per_transcript.tsv'))

        transcript_counts = defaultdict(int)

        bam = pysam.AlignmentFile(bam_file, 'rb')

        for read in bam.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            transcript = read.reference_name

            # if transcript.startswith('Chr_'):
            #     transcript = transcript[4:]

            transcript_counts[transcript] += 1

        bam.close()

        # --- Step 3: Write output ---
        with open(output_file, 'w') as out:

            out.write("Isoform_id\tGene_id\tRead_counts\n")

            for transcript, count in sorted(transcript_counts.items()):

                gene = transcript_to_gene.get(transcript, 'NA')
                
                out.write(f"{transcript}\t{gene}\t{count}\n")

        print(f'Saved: {output_file}')
                
# ============================================================
# COMPLETE WORKFLOW
# ============================================================

def run_complete_workflow(args):

    bambu_read_assignment(
        args.bam_files,
        args.output_dirs,
        args.genome_fasta,
        args.gtf
    )

    sort_fastq_reads_per_isoform(
        args.fastq_files,
        args.output_dirs
    )

    align_reads_to_isoform(
        args.output_dirs,
        args.bam_name,
        args.transcript_dir,
        args.threads,
        args.sam_flag,
        args.delete_intermediate
    )

    count_reads_per_transcript(
        args.output_dirs,
        args.bam_name,
        args.countby,
        args.gtf
    )


# ============================================================
# COMMAND LINE INTERFACE
# ============================================================

def parse_arguments():

    parser = argparse.ArgumentParser(
        description='Isoform-specific analysis pipeline'
    )

    subparsers = parser.add_subparsers(dest='command', required=True)

    prepare = subparsers.add_parser('prepare')
    prepare.add_argument('--gtf', required=True)
    prepare.add_argument('--genome_fasta', required=True)
    prepare.add_argument('--transcript_dir', required=True)

    assign = subparsers.add_parser('assign')
    assign.add_argument('--bam_files', nargs='+', required=True)
    assign.add_argument('--output_dirs', nargs='+', required=True)
    assign.add_argument('--genome_fasta', required=True)
    assign.add_argument('--gtf', required=True)

    sortfq = subparsers.add_parser('sort-fastq')
    sortfq.add_argument('--fastq_files', nargs='+', required=True)
    sortfq.add_argument('--output_dirs', nargs='+', required=True)

    align = subparsers.add_parser('align')
    align.add_argument('--output_dirs', nargs='+', required=True)
    align.add_argument('--transcript_dir', required=True)
    align.add_argument('--bam_name', default='filt_q8.bam')
    align.add_argument('--threads', type=int, default=4)
    align.add_argument('--sam_flag', type=int, default=2308)
    align.add_argument('--delete_intermediate', action='store_true')

    count = subparsers.add_parser('count')
    count.add_argument('--output_dirs', nargs='+', required=True)
    count.add_argument('--bam_name', required=True)
    count.add_argument('--gtf', required=True)
    count.add_argument('--countby', default='gene_name')

    all_cmd = subparsers.add_parser('all')
    all_cmd.add_argument('--bam_files', nargs='+', required=True)
    all_cmd.add_argument('--fastq_files', nargs='+', required=True)
    all_cmd.add_argument('--output_dirs', nargs='+', required=True)
    all_cmd.add_argument('--genome_fasta', required=True)
    all_cmd.add_argument('--gtf', required=True)
    all_cmd.add_argument('--transcript_dir', required=True)
    all_cmd.add_argument('--bam_name', default='filt_q8.bam')
    all_cmd.add_argument('--countby', default='gene_name')
    all_cmd.add_argument('--threads', type=int, default=4)
    all_cmd.add_argument('--sam_flag', type=int, default=2308)
    all_cmd.add_argument('--delete_intermediate', action='store_true')

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_arguments()

    if args.command == 'prepare':
        prepare_transcriptome(
            args.gtf,
            args.genome_fasta,
            args.transcript_dir
        )

    elif args.command == 'assign':
        bambu_read_assignment(
            args.bam_files,
            args.output_dirs,
            args.genome_fasta,
            args.gtf
        )

    elif args.command == 'sort-fastq':
        sort_fastq_reads_per_isoform(
            args.fastq_files,
            args.output_dirs
        )

    elif args.command == 'align':
        align_reads_to_isoform(
            args.output_dirs,
            args.bam_name,
            args.transcript_dir,
            args.threads,
            args.sam_flag,
            args.delete_intermediate
        )

    elif args.command == 'count':
        count_reads_per_transcript(
            args.output_dirs,
            args.bam_name,
            args.countby,
            args.gtf
        )

    elif args.command == 'all':
        run_complete_workflow(args)


if __name__ == '__main__':
    main()
