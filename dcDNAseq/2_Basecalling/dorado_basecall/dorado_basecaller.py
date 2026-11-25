#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dorado basecalling and post-processing pipeline
Refactored for reproducibility and reuse
"""

import os
import shutil
import argparse
import subprocess
from time import sleep


# ---------------- Utility ---------------- #

def run_command(cmd):
    print(f"Running: {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)


def extract_number(filename):
    base = os.path.basename(filename)
    num_part = base.split('_')[-1].split('.')[0]
    return int(num_part)


# ---------------- Dorado Basecalling ---------------- #

def run_dorado(input_dir, output_file, model, status_file, output_dir):

    status_path = os.path.join(output_dir, status_file)

    if not os.path.exists(status_path):

        cmd = f"dorado basecaller {model} {input_dir} > {output_file}"
        run_command(cmd)

        # Mark completion
        with open(status_path, "w") as f:
            f.write("done")

    else:
        print(f"{status_file} exists. Skipping.", flush=True)


# ---------------- Batch Processing ---------------- #

def run_dorado_batch(source_dirs, tmp_dir, output_dir,
                    model, file_type, batch_size):

    os.makedirs(output_dir, exist_ok=True)

    for source_dir in source_dirs:

        files = sorted(
            [f for f in os.listdir(source_dir) if f.endswith(file_type)],
            key=extract_number
        )

        batch_num = 1

        for i in range(0, len(files), batch_size):

            group = files[i:i + batch_size]

            bam_file = os.path.join(
                output_dir,
                f"{os.path.basename(source_dir)}_{batch_num}.bam"
            )

            status_file = f"{os.path.basename(source_dir)}_{batch_num}_done"

            # Clean tmp dir
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)

            if not os.path.exists(os.path.join(output_dir, status_file)):
                os.makedirs(tmp_dir, exist_ok=True)

                for f in group:
                    shutil.copy(
                        os.path.join(source_dir, f),
                        tmp_dir
                    )

            print(f"Processing batch {batch_num}: {group}")

            run_dorado(tmp_dir, bam_file, model, status_file, output_dir)

            batch_num += 1
            sleep(1)

        print(f"Finished processing {source_dir}")


# ---------------- Post-processing ---------------- #

def merge_bam(output_dir):

    files = os.listdir(output_dir)
    bam_files = sorted(
        [os.path.join(output_dir, f) for f in files if f.endswith(".bam")]
    )

    merged_dir = os.path.join(output_dir, "merged")
    os.makedirs(merged_dir, exist_ok=True)

    out_bam = os.path.join(merged_dir, "All_merged.bam")

    cmd = f"samtools merge -o {out_bam} " + " ".join(bam_files)
    run_command(cmd)


def get_summary_file(output_dir):

    bam = os.path.join(output_dir, "merged", "All_merged.bam")
    tsv = os.path.join(output_dir, "merged", "summary.tsv")

    cmd = f"dorado summary {bam} > {tsv}"
    run_command(cmd)


def bam_to_fastq(output_dir):

    bam = os.path.join(output_dir, "merged", "All_merged.bam")
    fastq = os.path.join(output_dir, "merged", "All_merged.fastq.gz")

    cmd = f"samtools bam2fq {bam} | pigz > {fastq}"
    run_command(cmd)


# ---------------- CLI ---------------- #

def main():

    parser = argparse.ArgumentParser(
        description="Dorado basecalling pipeline"
    )

    sub = parser.add_subparsers(dest="command")

    run_parser = sub.add_parser("run")
    run_parser.add_argument("--input_dirs", nargs="+", required=True)
    run_parser.add_argument("--tmp_dir", required=True)
    run_parser.add_argument("--output_dir", required=True)
    run_parser.add_argument("--model", required=True)
    run_parser.add_argument("--file_type", default=".pod5")
    run_parser.add_argument("--batch_size", type=int, default=1)

    args = parser.parse_args()

    if args.command == "run":

        run_dorado_batch(
            args.input_dirs,
            args.tmp_dir,
            args.output_dir,
            args.model,
            args.file_type,
            args.batch_size
        )

        merge_bam(args.output_dir)
        get_summary_file(args.output_dir)
        bam_to_fastq(args.output_dir)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()