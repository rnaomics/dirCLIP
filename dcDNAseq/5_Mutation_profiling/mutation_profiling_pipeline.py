#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mutation Profiling Pipeline

Subcommands:
    perbase
    mutation
    all
"""

import argparse
import gzip
import os
import subprocess
from pathlib import Path
import pandas as pd


def create_folder(directory):
    Path(directory).mkdir(parents=True, exist_ok=True)

def create_result_subdirs(result_dirs, min_read_quality, min_base_quality):
    """
    Create quality-specific output directories.

    Example:
    rep1_1pct/ -> rep1_1pct/q_8_Q_8/
    """

    return [
        os.path.join(
            result_dir,
            f"q_{min_read_quality}_Q_{min_base_quality}"
        )
        for result_dir in result_dirs
    ]

def calculate_mutation(mut_txt, mutation_tsv, inc_del, cov_thresh, mut_thresh):
    with gzip.open(mut_txt, 'rt') as infile, open(mutation_tsv, 'w') as outfile:
        for line in infile:
            if line.startswith('REF'):
                continue

            fields = line.strip().split()
            REF, POS, REF_BASE = fields[0], fields[1], fields[2]
            A, C, G, T = map(int, fields[4:8])
            DEL = int(fields[10])

            acgt = [A, C, G, T]
            nuc_count = sum(acgt)

            ref_idx = [REF_BASE.upper() == b for b in ['A','C','G','T']]
            ref_counts = sum(v for v, keep in zip(acgt, ref_idx) if keep)

            mutations = nuc_count - ref_counts
            coverage = nuc_count

            if inc_del:
                coverage += DEL
                mutations += DEL

            if coverage >= cov_thresh and mutations >= mut_thresh:
                outfile.write(
                    f"{REF}\t{POS}\t{coverage}\t{mutations/coverage}\n"
                )

    subprocess.call(f'gzip -f {mutation_tsv}', shell=True)


def get_filenames(min_read_quality, cov_thresh, mut_thresh):
    bam_name = f'filt_q{min_read_quality}.bam' if min_read_quality > 0 else 'sorted.bam'

    return {
        'filt_bam_file_name': bam_name,
        'perbase_output_file': 'mutation.txt',
        'coverage_tsv': f'Coverage_{cov_thresh}.tsv',
        'mutation_tsv': f'Mutation_rate_cov_{cov_thresh}_mut_{mut_thresh}.tsv'
    }


def run_perbase(bam_dirs, result_dirs, genome, filenames, max_depth, min_base_quality, threads):
    for bam_dir, result_dir in zip(bam_dirs, result_dirs):
        create_folder(result_dir)

        bam_file = os.path.join(bam_dir, filenames['filt_bam_file_name'])
        perbase_file = os.path.join(result_dir, filenames['perbase_output_file'])

        cmd = (
            f'perbase base-depth -t {threads} -D {max_depth} '
            f'-Q {min_base_quality} -r {genome} {bam_file} > {perbase_file}'
        )

        print(cmd)
        subprocess.call(cmd, shell=True)
        subprocess.call(f'gzip -f {perbase_file}', shell=True)


def run_calculate_mutation(result_dirs, filenames, inc_del, cov_thresh, mut_thresh):
    
    subdir = '2_Inc_del' if inc_del else '1_No_del'

    for result_dir in result_dirs:
        outdir = os.path.join(result_dir, subdir)
        create_folder(outdir)

        perbase_file = os.path.join(result_dir, filenames['perbase_output_file'] + '.gz')
        mutation_file = os.path.join(outdir, filenames['mutation_tsv'])

        calculate_mutation(perbase_file, mutation_file, inc_del, cov_thresh, mut_thresh)


def run_complete_workflow(args):
    filenames = get_filenames(args.min_read_quality, args.cov_thresh, args.mut_thresh)

    result_dirs = create_result_subdirs(args.output_dirs, args.min_read_quality, args.min_base_quality)

    run_perbase(args.bam_dirs, result_dirs, args.genome, filenames,
                args.max_depth, args.min_base_quality, args.threads)

    run_calculate_mutation(result_dirs, filenames,
                           args.incdel,
                           args.cov_thresh,
                           args.mut_thresh)


def build_parser():
    parser = argparse.ArgumentParser(description='Mutation profiling pipeline')
    sub = parser.add_subparsers(dest='command', required=True)

    def add_common(p):
        p.add_argument('--cov_thresh', type=int, default=10)
        p.add_argument('--mut_thresh', type=int, default=0)
        p.add_argument('--min_read_quality', type=int, default=8)
        p.add_argument('--min_base_quality', type=int, default=8)

    p1 = sub.add_parser('perbase')
    p1.add_argument('--bam_dirs', nargs='+', required=True)
    p1.add_argument('--output_dirs', nargs='+', required=True)
    p1.add_argument('--genome', required=True)
    p1.add_argument('--max_depth', type=int, default=500000)
    p1.add_argument('--threads', type=int, default=2)
    add_common(p1)

    p3 = sub.add_parser('mutation')
    p3.add_argument('--output_dirs', nargs='+', required=True)
    p3.add_argument('--incdel', action='store_true', help='Include deletions in mutation rate calculation')
    add_common(p3)

    p4 = sub.add_parser('all')
    p4.add_argument('--bam_dirs', nargs='+', required=True)
    p4.add_argument('--output_dirs', nargs='+', required=True)
    p4.add_argument('--genome', required=True)
    p4.add_argument('--max_depth', type=int, default=500000)
    p4.add_argument('--threads', type=int, default=2)
    p4.add_argument('--incdel', action='store_true', help='Include deletions in mutation rate calculation')
    add_common(p4)

    return parser


def main():
    args = build_parser().parse_args()

    if args.command == 'all':
        run_complete_workflow(args)

    elif args.command == 'perbase':
        filenames = get_filenames(args.min_read_quality, args.cov_thresh, args.mut_thresh)
        result_dirs = create_result_subdirs(args.output_dirs, args.min_read_quality, args.min_base_quality)
        run_perbase(args.bam_dirs, result_dirs, args.genome, filenames,
                    args.max_depth, args.min_base_quality, args.threads)

    elif args.command == 'mutation':
        filenames = get_filenames(args.min_read_quality, args.cov_thresh, args.mut_thresh)
        result_dirs = create_result_subdirs(args.output_dirs, args.min_read_quality, args.min_base_quality)
        run_calculate_mutation(result_dirs, filenames,
                               args.incdel,
                               args.cov_thresh,
                               args.mut_thresh)


if __name__ == '__main__':
    main()
