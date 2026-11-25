#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Minimap2 alignment and BAM post-processing pipeline
"""

import os
import argparse
import subprocess
import pysam


# ---------------- Utilities ---------------- #

def run_command(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def create_folder(path):
    os.makedirs(path, exist_ok=True)


# ---------------- Alignment ---------------- #

def run_minimap2(fastq_dirs, genome, ss_file,
                 output_dirs, bam_name, threads):

    for fq_dir, out_dir in zip(fastq_dirs, output_dirs):

        create_folder(out_dir)

        fastq_path = fq_dir
        bam_path = os.path.join(out_dir, bam_name)

        
        # --- Conditional minimap2 command ---
        if ss_file:
            mm2_cmd = (
                f"minimap2 -ax splice -uf -t {threads} -k14 "
                f"--junc-bed {ss_file} {genome} {fastq_path}"
            )
        else:
            mm2_cmd = (
                f"minimap2 -ax splice -uf -t {threads} -k14 "
                f"{genome} {fastq_path}"
            )

        # --- Full pipeline command ---
        cmd = (
            f"{mm2_cmd} 2> {out_dir}/error.log | "
            f"samtools view -bh -F 2308 | "
            f"samtools sort -O BAM -o {bam_path}"
        )

        run_command(cmd)

        bam_index(out_dir, bam_name)
        alignment_stats(out_dir, bam_name)
        alignment_flags(out_dir, bam_name)
        bam_to_bed(out_dir, bam_name)


# ---------------- BAM Processing ---------------- #

def bam_index(out_dir, bam_name):
    run_command(f"samtools index {os.path.join(out_dir, bam_name)}")


def alignment_stats(out_dir, bam_name):
    stats_file = bam_name.replace(".bam", "_alignmentStats.txt")
    run_command(
        f"samtools flagstat {os.path.join(out_dir, bam_name)} "
        f"> {os.path.join(out_dir, stats_file)}"
    )


def alignment_flags(out_dir, bam_name):
    stats_file = bam_name.replace(".bam", "_flag_stats.txt")

    run_command(
        f"samtools view {os.path.join(out_dir, bam_name)} "
        f"| awk '{{print $2}}' > {os.path.join(out_dir, stats_file)}"
    )


def bam_to_bed(out_dir, bam_name):
    bed_file = bam_name.replace(".bam", ".bed")
    bed_path = os.path.join(out_dir, bed_file)

    run_command(f"bedtools bamtobed -i {os.path.join(out_dir, bam_name)} > {bed_path}")
    run_command(f"gzip -f {bed_path}")


# ---------------- Filtering ---------------- #

def filter_bam_quality(output_dirs, min_q, bam_name):

    if min_q <= 0:
        return

    filt_name = f"filt_q{min_q}.bam"

    for out_dir in output_dirs:
        input_bam = os.path.join(out_dir, bam_name)
        output_bam = os.path.join(out_dir, filt_name)

        with pysam.AlignmentFile(input_bam, "rb") as infile, \
             pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:

            for read in infile:
                if read.query_qualities:
                    mean_q = sum(read.query_qualities) / len(read.query_qualities)
                    if mean_q > min_q:
                        outfile.write(read)

        bam_index(out_dir, filt_name)
        alignment_stats(out_dir, filt_name)
        bam_to_bed(out_dir, filt_name)


# ---------------- CLI ---------------- #

def main():
    parser = argparse.ArgumentParser(description="Minimap2 alignment pipeline")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # -------- ALIGN subcommand -------- #
    align_parser = subparsers.add_parser("align", help="Run minimap2 alignment")

    align_parser.add_argument("--fastq_dirs", nargs="+", required=True)
    align_parser.add_argument("--output_dirs", nargs="+", required=True)

    align_parser.add_argument("--genome", required=True)
    align_parser.add_argument("--splice_sites", default=None, help="Splice junction BED file (optional)")

    align_parser.add_argument("--bam_name", default="sorted.bam")
    align_parser.add_argument("--threads", type=int, default=4)

    # -------- FILTER subcommand -------- #
    filter_parser = subparsers.add_parser("filter", help="Filter BAM by quality")

    filter_parser.add_argument("--output_dirs", nargs="+", required=True)
    filter_parser.add_argument("--bam_name", default="sorted.bam")
    filter_parser.add_argument("--min_quality", type=int, default=0)

    args = parser.parse_args()

    # -------- Execution -------- #
    if args.command == "align":
        run_minimap2(
            args.fastq_dirs,
            args.genome,
            args.splice_sites,
            args.output_dirs,
            args.bam_name,
            args.threads
        )

    elif args.command == "filter":
        filter_bam_quality(
            args.output_dirs,
            args.min_quality,
            args.bam_name
        )




if __name__ == "__main__":
    main()