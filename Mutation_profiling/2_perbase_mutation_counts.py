# -*- coding: utf-8 -*-
"""
@author: mohamedb
"""

import pandas as pd
import subprocess
import os
import gzip
import pysam

################ 1. Create folder if it doesn't exist ################
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

################ 2. Calculate mutations from mutation file ################
def calculate_mutation(mut_txt, mut_bed, inc_del, inc_ins, cov_thresh, Mut_thresh):
    with gzip.open(mut_txt, 'rt') as wiggle_file, open(mut_bed, 'w') as bed_file:
        for line in wiggle_file:
            if line.startswith('REF'):
                continue

            REF, POS, REF_BASE, DEPTH, A, C, G, T, N, INS, DEL, REF_SKIP, FAIL, NEAR_MAX_DEPTH = line.strip().split()
            ACGT = [int(A), int(C), int(G), int(T)]
            Nuc_count = sum(ACGT)

            id = [REF_BASE.upper() == i for i in ['A', 'C', 'G', 'T']]
            Ref_counts = sum([value for value, include in zip(ACGT, id) if include])
            Mutations = Nuc_count - Ref_counts
            Coverage = Nuc_count

            if inc_ins:
                Mutations += int(INS)
            if inc_del:
                Coverage += int(DEL)
                Mutations += int(DEL)

            if Coverage >= cov_thresh and Mutations >= Mut_thresh:
                Mutation_rate = Mutations / Coverage
                bed_file.write(f'{REF}\t{POS}\t{Coverage}\t{Mutation_rate}\n')

    subprocess.call(f"gzip -f {mut_bed}", shell=True)

################ 3. Run perbase tool on BAM files ################
def run_perbase(Bam_dir, Result_dir, genomes, filenames, max_depth, min_base_quality_score):
    bam_file_name = filenames["filt_bam_file_name"]
    mut_file_name = filenames["mut_file_name"]

    for i in range(len(Bam_dir)):
        genome = genomes[i]
        genome_file = {
            "mm9": Project_dir + "genome/mm9/mm9.fa",
            "mm39": Project_dir + "genome/mm39/mm39.fa",
            "7sk": Project_dir + "genome/7SK_RNA/7sk_reverse.fasta"
        }.get(genome)

        bam_file = Bam_dir[i] + bam_file_name
        mut_file = Result_dir[i] + mut_file_name
        createFolder(Result_dir[i])

        command = f"perbase base-depth -t 2 -D {max_depth} -Q {min_base_quality_score} -r {genome_file} {bam_file} > {mut_file}"
        print(command)
        subprocess.call(command, shell=True)

        subprocess.call(f"gzip -f {mut_file}", shell=True)

################ 4. Run mutation calculation ################
def run_calculate_mutation(Result_dir, filenames, mut_flags, cov_thresh, mut_thresh):
    sub_res_dir = mut_flags["sub_res_dir"]
    Inc_del = mut_flags["Inc_del"]
    Inc_ins = mut_flags["Inc_ins"]

    mut_file_name = filenames["mut_file_name"]
    mut_bed_file_name = filenames["mut_bed_file_name"]

    for i in range(len(Result_dir)):
        res_path = Result_dir[i] + sub_res_dir
        createFolder(res_path)

        mut_file = Result_dir[i] + mut_file_name + ".gz"
        mut_bed = res_path + mut_bed_file_name

        calculate_mutation(mut_file, mut_bed, Inc_del, Inc_ins, cov_thresh, mut_thresh)

################ 5. Calculate per-nucleotide coverage ################
def calculate_coverage(mut_txt, cov_bed, cov_thresh):
    if os.path.exists(cov_bed):
        subprocess.call(f"rm {cov_bed}", shell=True)

    with open(cov_bed, 'a') as bed_file:
        for chunk in pd.read_csv(mut_txt, compression='gzip', sep='\t', chunksize=1_000_000):
            chunk['Coverage'] = chunk[['A', 'C', 'G', 'T', 'N', 'DEL']].sum(axis=1)
            chunk_filtered = chunk[chunk['Coverage'] >= cov_thresh][['REF', 'POS', 'Coverage']]
            chunk_filtered.to_csv(bed_file, sep='\t', index=False, header=False)

    subprocess.call(f"gzip -f {cov_bed}", shell=True)

################ 6. Run coverage calculation ################
def run_calculate_coverage(Result_dir, filenames, cov_thresh):
    mut_file_name = filenames["mut_file_name"]
    cov_bed_file_name = filenames["cov_bed_file_name"]

    for i in range(len(Result_dir)):
        res_path = Result_dir[i]
        mut_file = Result_dir[i] + mut_file_name + ".gz"
        cov_bed = res_path + cov_bed_file_name
        createFolder(res_path)
        calculate_coverage(mut_file, cov_bed, cov_thresh)

################ 7. Filter BAM files by read quality ################
def filter_bam_by_read_quality(Bam_dir, min_read_quality, bam_file_name):
    if min_read_quality > 0:
        filt_bam_file_name = "filt_q" + str(min_read_quality) + ".bam"
        for i in range(len(Bam_dir)):
            input_bam = Bam_dir[i] + bam_file_name
            output_bam = Bam_dir[i] + filt_bam_file_name

            with pysam.AlignmentFile(input_bam, "rb") as infile, \
                 pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
                for read in infile:
                    if read.query_qualities:
                        mean_quality = sum(read.query_qualities) / len(read.query_qualities)
                        if mean_quality > min_read_quality:
                            outfile.write(read)

            subprocess.call(f"samtools index {output_bam}", shell=True)

################ 8. Set mutation flags ################
def set_mut_flags(mut_criteria):
    flags = {
        "Inc_del": {"Inc_del": True, "Inc_ins": False, "sub_res_dir": "2_Inc_del/"},
        "No_del": {"Inc_del": False, "Inc_ins": False, "sub_res_dir": "1_No_del/"},
        "Inc_ins_del": {"Inc_del": True, "Inc_ins": True, "sub_res_dir": "3_Inc_ins_del/"}
    }
    return flags.get(mut_criteria)

################ 9. Build filenames ################
def get_filenames(min_read_quality, cov_thresh, mut_thresh):
    bam_file_name = "sorted.bam"
    filt_bam_file_name = f"filt_q{min_read_quality}.bam" if min_read_quality > 0 else bam_file_name
    mut_file_name = "mutation.txt"
    cov_bed_file_name = f"Coverage_{cov_thresh}.tsv"
    mut_bed_file_name = f"Mutation_rate_cov_{cov_thresh}_mut_{mut_thresh}.tsv"

    return {
        "bam_file_name": bam_file_name,
        "filt_bam_file_name": filt_bam_file_name,
        "mut_file_name": mut_file_name,
        "cov_bed_file_name": cov_bed_file_name,
        "mut_bed_file_name": mut_bed_file_name
    }

################ ###### Pipeline Execution Block ################
Project_dir = "/data/groups/anko/mohamed/oxford_nanopore/"

Bam_dir = ["Results/Nanopore/5_UV_transcriptome/Alignment/minimap/mm9/Mauro/TotalRNA_cDNA_UV_merge/"]
sub_dirs = [f"subsample{i}/" for i in range(1, 4)]
Result_dir = ["Results/Nanopore/5_UV_transcriptome/Mutation_profiling/minimap/mm9/Mauro/TotalRNA_cDNA_UV_merge/"]
genomes = ["mm9"]
cov_thresh = 10
mut_thresh = 0
min_base_quality_scores = [8]
min_read_qualitys = [8]
mut_criterias = ["No_del", "Inc_del"]

os.chdir(Project_dir)

# Combine subsamples with base directories
Bam_dirs = [b + s for b in Bam_dir for s in sub_dirs]
genomes2 = genomes * len(sub_dirs)

bam_filename = "sorted.bam"

for min_read_quality in set(min_read_qualitys):
    filter_bam_by_read_quality(Bam_dirs, min_read_quality, bam_filename)

for ii in range(len(min_base_quality_scores)):
    min_base_quality_score = min_base_quality_scores[ii]
    min_read_quality = min_read_qualitys[ii]
    filenames = get_filenames(min_read_quality, cov_thresh, mut_thresh)

    Result_dirs = [r + f"q_{min_read_quality}_Q_{min_base_quality_score}/" + s for r in Result_dir for s in sub_dirs]

    run_perbase(Bam_dirs, Result_dirs, genomes2, filenames, max_depth=500000, min_base_quality_score=min_base_quality_score)

    for mut_criteria in mut_criterias:
        mut_flags = set_mut_flags(mut