#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Nanopore basecalling and FASTQ merging pipeline

Author: Mohamed B.
Refactored for reproducibility and reuse.
"""

import os
import argparse
import subprocess


# ---------------- Utility Functions ---------------- #

def create_folder(directory):
    """Create directory if it does not exist."""
    os.makedirs(directory, exist_ok=True)


def run_command(cmd):
    """Run shell command with logging."""
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


# ---------------- Basecalling ---------------- #

def guppy_basecall(source_dirs, flowcell, kit, callers, threads):
    """Run Guppy basecalling on input directories."""
    
    for source_dir in source_dirs:
        fast5_dirs = [
            os.path.join(source_dir, "fast5_fail"),
            os.path.join(source_dir, "fast5_pass"),
            os.path.join(source_dir, "fast5_skip"),
        ]

        out_dirs = [
            os.path.join(source_dir, "guppy_out", "fail"),
            os.path.join(source_dir, "guppy_out", "pass"),
            os.path.join(source_dir, "guppy_out", "skip"),
        ]

        for fast5_dir, out_dir in zip(fast5_dirs, out_dirs):
            create_folder(out_dir)

            log_file = os.path.join(out_dir, "guppy.log")

            cmd = (
                f"guppy_basecaller "
                f"-i {fast5_dir} "
                f"-s {out_dir} "
                f"--flowcell {flowcell} "
                f"--kit {kit} "
                f"--num_callers {callers} "
                f"--cpu_threads_per_caller {threads} "
                f"-r > {log_file}"
            )

            run_command(cmd)


# ---------------- FASTQ Merging ---------------- #

def merge_fastq(source_dirs):
    """Merge FASTQ outputs from guppy."""

    for source_dir in source_dirs:

        dest_dir = os.path.join(source_dir, "merge_guppy_outputs")
        create_folder(dest_dir)

        pass_dirs = [
            os.path.join(source_dir, "guppy_out", x, "pass")
            for x in ["fail", "pass", "skip"]
        ]

        fail_dirs = [
            os.path.join(source_dir, "guppy_out", x, "fail")
            for x in ["fail", "pass", "skip"]
        ]

        pass_merge_files = [
            os.path.join(dest_dir, f"{x}_pass.fastq")
            for x in ["fail", "pass", "skip"]
        ]

        fail_merge_files = [
            os.path.join(dest_dir, f"{x}_fail.fastq")
            for x in ["fail", "pass", "skip"]
        ]

        # Merge PASS
        for d, out in zip(pass_dirs, pass_merge_files):
            run_command(f"cat {d}/*.fastq > {out}")

        run_command(
            f"cat {' '.join(pass_merge_files)} > {dest_dir}/All_pass.fastq"
        )

        # Merge FAIL
        for d, out in zip(fail_dirs, fail_merge_files):
            run_command(f"cat {d}/*.fastq > {out}")

        run_command(
            f"cat {' '.join(fail_merge_files)} > {dest_dir}/All_fail.fastq"
        )

        # Merge ALL
        run_command(
            f"cat {dest_dir}/All_fail.fastq {dest_dir}/All_pass.fastq "
            f"> {dest_dir}/All_merged.fastq"
        )

        # Clean intermediate files
        for f in pass_merge_files + fail_merge_files:
            if os.path.exists(f):
                os.remove(f)


# ---------------- CLI ---------------- #

def main():
    parser = argparse.ArgumentParser(
        description="Guppy basecalling"
    )

    subparsers = parser.add_subparsers(dest="command")

    # Basecall command
    basecall_parser = subparsers.add_parser("basecall")
    basecall_parser.add_argument("--input_dirs", nargs="+", required=True)
    basecall_parser.add_argument("--flowcell", default="FLO-MIN106")
    basecall_parser.add_argument("--kit", default="SQK-RNA002")
    basecall_parser.add_argument("--callers", type=int, default=2)
    basecall_parser.add_argument("--threads", type=int, default=5)

    # Merge command
    merge_parser = subparsers.add_parser("merge")
    merge_parser.add_argument("--input_dirs", nargs="+", required=True)

    args = parser.parse_args()

    if args.command == "basecall":
        guppy_basecall(
            args.input_dirs,
            args.flowcell,
            args.kit,
            args.callers,
            args.threads,
        )

    elif args.command == "merge":
        merge_fastq(args.input_dirs)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()