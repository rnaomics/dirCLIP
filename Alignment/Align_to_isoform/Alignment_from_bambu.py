# conda activate isoquant

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 2022
@author: mohamedb
"""

import os
import shutil
import subprocess
import glob
import pandas as pd
import pdb
import pysam
from collections import defaultdict

# 📁 Utility: Create directory if it doesn't exist
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

# 📦 Move source file to destination
def mov_files(source_file, dest_file):
    command = "mv -f " + source_file + " " + dest_file
    print(command)
    subprocess.call(command, shell=True)

# 🧬 Run Bambu isoform-specific read assignment
def bambu_read_assignment(Bam_files, bambu_dirs, genome_file, gtf_file, Project_dir, Tools_dir):
    os.chdir(Tools_dir)

    # Write parameters to file
    with open('parameters.txt', 'w') as f:
        f.write(f"Project_dir={Project_dir}\n")
        f.write(f"genome_file={genome_file}\n")
        f.write(f"gtf_file={gtf_file}\n")
        f.write(f"Bam_files={','.join(Bam_files)}\n")
        f.write(f"bambu_dirs={','.join(bambu_dirs)}\n")

    # Execute R script
    command = "R CMD BATCH bambu_for_read_assignment.R"
    subprocess.call(command, shell=True)

    # Clean up temporary files
    remove_files('.RData')
    # remove_files('parameters.txt')  # Optional

    os.chdir(Project_dir)

# 🧬 Add gene IDs to isoform counts
def add_geneid_to_isoform_counts(transcriptid_gene_id_path, isoform_counts_file):
    Transcript_gene_id = pd.read_csv(transcriptid_gene_id_path, sep='\t')
    Transcript_gene_id.set_index('transcript_id', inplace=True)

    isoform_counts = pd.read_csv(isoform_counts_file, header=None, sep=';')
    isoform_counts.columns = ["Isoform_id", "Read_counts"]

    isoform_counts['Gene_id'] = isoform_counts['Isoform_id'].map(Transcript_gene_id['gene_id'])
    isoform_counts['N_isoforms'] = isoform_counts.groupby('Gene_id')['Isoform_id'].transform('count')

    isoform_counts = isoform_counts[['Gene_id', 'N_isoforms', 'Isoform_id', 'Read_counts']]
    isoform_counts.sort_values(by='Gene_id', inplace=True)

    output_file = isoform_counts_file.replace(".csv", "_with_geneid.csv")
    isoform_counts.to_csv(output_file, sep=';', index=False)

############## 2. Function to sort fastq reads per isoform ###############################################
def sort_fastq_reads_per_isoform(sort_fastq, fastq_files, dirs, transcriptid_gene_id_path):
    for i in range(len(fastq_files)):
        fastq_infile = fastq_files[i]
        fastq_out = dirs[i] + "fastq_by_isoform/"
        read_assignment_file = dirs[i] + "Out/readToTranscriptMaps.tsv"

        createFolder(fastq_out)

        command = "python " + sort_fastq + " -fi " + fastq_infile + " -fo " + fastq_out + \
            " -ir " + read_assignment_file
        
        print(command)
        subprocess.call(command, shell=True)

        ### move the isoform_counts to another folder
        source_file = dirs[i] + "fastq_by_isoform/isoform_counts.csv"
        dest_file = dirs[i] + "isoform_counts.csv"

        command = "mv -f " + source_file + " " + dest_file
        print(command)
        subprocess.call(command, shell=True)

        #### Adding gene_id to the isoform_counts ###########
        add_geneid_to_isoform_counts(transcriptid_gene_id_path, dest_file)

# 🧹 Remove a file or folder
def remove_files(Path):
    command = "rm -r " + Path
    subprocess.call(command, shell=True)

# 🔀 Merge BAM files in batches
def merge_bam_files_in_batches(bam_directory, output_file, batch_size=100):
    bam_files = [os.path.join(bam_directory, f) for f in os.listdir(bam_directory) if f.endswith(".bam")]
    bam_files.sort()

    batches = [bam_files[i:i + batch_size] for i in range(0, len(bam_files), batch_size)]
    intermediate_files = []

    for i, batch in enumerate(batches):
        partial_output = output_file.replace("merged.bam", f"partial_merged_{i}.bam")
        intermediate_files.append(partial_output)

        command = ["samtools", "merge", "-f", partial_output] + batch
        print(f"Merging batch {i+1}/{len(batches)}: {command}")
        subprocess.call(command)

    final_merge_command = ["samtools", "merge", "-f", output_file] + intermediate_files
    print(f"Merging all partial files into {output_file}: {final_merge_command}")
    subprocess.call(final_merge_command)

    # Replace header
    filt_header = output_file.replace("merged.bam", "filt_header.txt")
    new_merged_bam = output_file.replace("merged.bam", "new_merged.bam")

    command = "samtools view -H " + output_file + " | grep -E '^@HD|^@SQ' > " + filt_header
    subprocess.call(command, shell=True)

    command = "samtools reheader " + filt_header + " " + output_file + " > " + new_merged_bam
    subprocess.call(command, shell=True)

    command = "mv " + new_merged_bam + " " + output_file
    subprocess.call(command, shell=True)

    os.remove(filt_header)

    for intermediate_file in intermediate_files:
        os.remove(intermediate_file)
        print(f"Removed intermediate file: {intermediate_file}")

# 🧬 Align reads to isoform-specific references
def align_reads_to_isoform(dirs, bam_file, fasta_path, Threads, sam_flag, delete_intermediate_files):
    for i in range(len(dirs)):
        fastq_path = dirs[i] + "fastq_by_isoform/"
        bam_path = dirs[i] + "bam_by_isoform/"

        createFolder(bam_path)

        isoform_counts_file = dirs[i] + "isoform_counts.csv"
        isoform_counts = pd.read_csv(isoform_counts_file, header=None)
        isoform_counts[['isoform', 'counts']] = isoform_counts[0].str.split(';', expand=True)
        isoform_counts = isoform_counts.drop(columns=[0])
        isoform_counts['counts'] = pd.to_numeric(isoform_counts['counts'])

        for isoform in isoform_counts["isoform"]:
            fastq_file = fastq_path + isoform + ".fastq.gz"
            fasta_file = fasta_path + "Chr_" + isoform + ".fa"
            bamfile = bam_path + isoform + ".bam"

            command = (
                f"minimap2 -ax splice -uf -t {Threads} -k14 {fasta_file} {fastq_file} | "
                f"samtools view -bh -F {sam_flag} | "
                f"samtools sort -O BAM -o {bamfile}"
            )
            subprocess.call(command, shell=True)

        merged_bam = dirs[i] + "merged.bam"
        merge_bam_files_in_batches(bam_path, merged_bam, 1000)

        sorted_bam = dirs[i] + bam_file
        subprocess.call(f"samtools sort -o {sorted_bam} {merged_bam}", shell=True)
        subprocess.call(f"samtools index {sorted_bam}", shell=True)

        align_stats = dirs[i] + "alignmentStats.txt"
        subprocess.call(f"samtools flagstat {sorted_bam} >> {align_stats}", shell=True)

        bed_file = dirs[i] + "sorted.bed.gz"
        subprocess.call(f"bedtools bamtobed -i {sorted_bam} | pigz > {bed_file}", shell=True)

        if delete_intermediate_files:
            remove_files(fastq_path)
            remove_files(bam_path)
            remove_files(merged_bam)

############################### 🧮 Function to count reads from BAM files per transcript ###############################

def Count_bam_readcounts(bambu_dirs, bambu_bam, gene_id, gtf_file):
    # 📜 Step 1: Parse GTF to map transcript IDs to gene names
    transcript_to_gene = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] != "transcript":
                continue
            attributes = fields[8]
            attrs = {
                kv.strip().split(' ')[0]: kv.strip().split(' ')[1].strip('"')
                for kv in attributes.strip(';').split('; ') if ' ' in kv
            }
            transcript_id = attrs.get("transcript_id")
            gene_name = attrs.get(gene_id, "NA")
            if transcript_id:
                transcript_to_gene[transcript_id] = gene_name

    # 📊 Step 2: Count reads per transcript from each BAM file
    for bambu_dir in bambu_dirs:
        bam_file = bambu_dir + bambu_bam
        output_file = bam_file.replace(".bam", "_read_counts_per_transcript.tsv")

        bam = pysam.AlignmentFile(bam_file, "rb")
        transcript_counts = defaultdict(int)

        for read in bam.fetch(until_eof=True):
            if not read.is_unmapped:
                transcript = read.reference_name
                if transcript.startswith("Chr_"):
                    transcript = transcript.replace("Chr_", "", 1)
                transcript_counts[transcript] += 1

        bam.close()

        # 📄 Step 3: Write counts to output file
        with open(output_file, "w") as f:
            f.write("Isoform_id\tGene_id\tRead_counts\n")
            for transcript, count in sorted(transcript_counts.items()):
                gene = transcript_to_gene.get(transcript, "NA")
                f.write(f"{transcript}\t{gene}\t{count}\n")

        print(f"✅ Saved read counts with gene names to {output_file}")

###########################################################################################################
Tools_dir = os.getcwd()
Project_dir = "/data/groups/anko/mohamed/oxford_nanopore/"
os.chdir(Project_dir)  # Switch to project directory

####################################################################################################################################

# BAM and FASTQ file paths
Bam_files = [
    "Results/Nanopore/3_SRSF3_DIRCLIP/Alignment/minimap/mm9/Mauro/CLIP_cDNA_rep1/filt_q8.bam",
    "Results/Nanopore/3_SRSF3_DIRCLIP/Alignment/minimap/mm9/Mauro/CLIP_cDNA_rep2/filt_q8.bam",
    "Results/Nanopore/4_noUV_transcriptome/Alignment/minimap/mm9/Mauro/TotalRNA_cDNA_no_UV_rep3/filt_q8.bam"
]

fastq_files = [
    "Data/Nanopore/3_SRSF3_DIRCLIP/Mauro/CLIP_cDNA_rep1/merge_dorado_outputs/All_merged.fastq.gz",
    "Data/Nanopore/3_SRSF3_DIRCLIP/Mauro/CLIP_cDNA_rep2/merge_dorado_outputs/All_merged.fastq.gz",
    "Data/Nanopore/4_noUV_transcriptome/Mauro/TotalRNA_cDNA_no_UV_rep3/merge_dorado_outputs/All_merged.fastq.gz"
]

bambu_dirs = [
    "Results/Nanopore/3_SRSF3_DIRCLIP/Alignment/bambu/mm9_transcriptome/CLIP_cDNA_rep1/",
    "Results/Nanopore/3_SRSF3_DIRCLIP/Alignment/bambu/mm9_transcriptome/CLIP_cDNA_rep2/",
    "Results/Nanopore/4_noUV_transcriptome/Alignment/bambu/mm9_transcriptome/Mauro/TotalRNA_cDNA_no_UV_rep3/"
]

bambu_bam = "filt_q8.bam"

################ Genome, annotation and Transcripts #####################################################################################

genome_file = "genome/mm9/mm9.fa"

# RefSeq annotation (default)
gtf_file = "genome/mm9/Annotation/mm9.refGene.gtf"
Transcript_fasta_path = "Data/Transcriptome/mm9/Transcriptome_from_genome/Transcripts/"
transcriptid_gene_id_path = "Data/Transcriptome/mm9/Transcriptome_from_genome/mm9_transcriptid_to_geneid.tsv.gz"
gene_id = "gene_id"

# Gencode annotation (alternative — toggle by uncommenting block below)
"""
gtf_file = "genome/mm9/Annotation/gencode.vM1.annotation.gtf"
Transcript_fasta_path = "Data/Transcriptome/mm9/Transcriptome_from_genome/gencode/Transcripts/"
transcriptid_gene_id_path = "Data/Transcriptome/mm9/Transcriptome_from_genome/gencode/mm9_transcriptid_to_geneid.tsv.gz"
gene_id = "gene_name"
"""

######################## 0. Preprocessing steps splitting single transcriptome file into multiple transcriptome file ################

Path = "Tools/get_transcriptome_from_genome/2_get_mm9_transcriptome/"

# a. Add "Chr_" prefix to transcript IDs
subprocess.call(f"python {Path}1b_add_chr_to_transcriptid.py", shell=True)

# b. Map transcript IDs to gene IDs
subprocess.call(f"python {Path}2_get_transcriptid_to_geneid.py", shell=True)

# c. Extract transcriptome FASTAs from genome
subprocess.call(f"{Path}2_get_transcriptome_from_genome.sh", shell=True)

# d. Split transcriptome into per-transcript FASTA files
subprocess.call(f"python {Path}3_splitting_transcriptome_transcriptwise_fasta.py", shell=True)

######################## 1. Performing bambu isoform specific read assignment ####################################################
bambu_read_assignment(Bam_files, bambu_dirs, genome_file, gtf_file, Project_dir, Tools_dir)

######################## 2. Sorting the fastq reads based isoquant read assignment file #############################################
sort_fastq = Tools_dir + "/sort_fastq_Bambu.py"

sort_fastq_reads_per_isoform(sort_fastq, fastq_files, bambu_dirs, transcriptid_gene_id_path)

######################## 3. Performing isoform specific alignment ###################################################################
sam_flag = 2308
delete_intermediate_files = True 
Threads = 2

align_reads_to_isoform(bambu_dirs, bambu_bam, Transcript_fasta_path, Threads, sam_flag, delete_intermediate_files)

######################## 4. Count the number of reads per isoform #######################################################################
Count_bam_readcounts(bambu_dirs, bambu_bam, gene_id, gtf_file)

############################################################################################################################################