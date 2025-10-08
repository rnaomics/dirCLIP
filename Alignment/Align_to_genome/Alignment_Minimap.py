"""
🧬 Refactored Pipeline: Nanopore Transcriptomic Analysis

Author: Mohamed B Tools Used: fastQC, minimap2, samtools, bedtools, Rsubread (via Rscript) 
Purpose: Full pipeline for long-read RNA sequencing analysis using dRNA samples, with quality control, alignment, and gene-wise read counting

"""

import os
import shutil
import subprocess
import glob
from pickle import TRUE
from re import T
from ssl import get_protocol_name

# 🗂️ Create a directory if it doesn't exist
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

# 🧪 Perform FastQC on all FASTQ files in a directory
def fastqc_check(fastq_dir, Result_dir):
    createFolder(Result_dir)

    files = [os.path.basename(x) for x in glob.glob(fastq_dir + "*.fastq.gz")]
    print(files)

    for file in files:
        qc_command = f"fastqc -o {Result_dir} --extract -f fastq {fastq_dir}{file}"
        print(qc_command + "\n")
        subprocess.call(qc_command, shell=True)
    
    multiqc_command = f"multiqc {Result_dir}"
    subprocess.call(multiqc_command, shell=True)

# ⚡ Run FastQC on selected FASTQ across multiple directories using Singularity
def run_fastqc(fastq_dirs, fastq_select, Result_dirs, threads):
    fastqc_exe = "singularity exec /data/containers/biobase.sif fastqc"
    for i in range(len(fastq_dirs)):
        fastq_dir = fastq_dirs[i]
        res_path = Result_dirs[i]
        createFolder(res_path)

        qc_command = f"{fastqc_exe} -o {res_path} --extract -f fastq -t {threads} {fastq_dir}{fastq_select}"
        print(qc_command)
        subprocess.call(qc_command, shell=True)

# 🧬 Run minimap2 alignment and generate BAM, BED, and QC files
def run_minimap2_aligner(fastq_dir, fastq_select, genomes, Result_dir, bam_file, bed_file):
    for ii in range(len(genomes)):
        genome = genomes[ii]

        if genome == "mm9":
            genome_file = "genome/mm9/mm9.fa"
            Sub_path_align = "Alignment/minimap/mm9/"
            ss_file = "genome/mm9/Annotation/mm9.knownGene.ss.bed"
        elif genome == "mm39":
            genome_file = "genome/mm39/mm39.fa"
            Sub_path_align = "Alignment/minimap/mm39/"
            ss_file = "genome/mm9/Annotation/mm39.ss.bed"
        elif genome == "hg38":
            genome_file = "genome/hg38/hg38.fa"
            Sub_path_align = "Alignment/minimap/hg38/"
            ss_file = "genome/hg38/Annotation/hg38.ss.bed"

        Result_dir2 = [x.replace("Sub_path/", Sub_path_align) for x in Result_dir]
        samtools = "singularity exec /data/containers/biobase.sif samtools"

        for i in range(len(fastq_dir)):
            file_path = fastq_dir[i] + fastq_select
            res_path = Result_dir2[i]
            createFolder(res_path)

            align_cmd = (
                f"minimap2 -ax splice -uf -t 5 -k14 --junc-bed {ss_file} {genome_file} {file_path} 2> {res_path}error.log | "
                f"{samtools} view -bh -F 2308 | {samtools} sort -O BAM -o {res_path}{bam_file}"
            )
            print(align_cmd)
            subprocess.call(align_cmd, shell=True)

            subprocess.call(f"{samtools} index {res_path}{bam_file}", shell=TRUE)
            subprocess.call(f"{samtools} flagstat {res_path}{bam_file} >> {res_path}alignmentStats.txt", shell=TRUE)

            flag_cmd = f"{samtools} view {res_path}{bam_file} | awk '{{print $2}}' > {res_path}flag_stats.txt"
            print(flag_cmd)
            subprocess.call(flag_cmd, shell=TRUE)

            bed_cmd = f"bedtools bamtobed -i {res_path}{bam_file} > {res_path}{bed_file}"
            print(bed_cmd)
            subprocess.call(bed_cmd, shell=True)

# 🔤 Convert integers to strings (for parameters)
def int_2_str(input):
    return [str(x) for x in input] if isinstance(input, list) else str(input)

# 🧾 Write parameters to file
def write_parameters(Parameters):
    with open('parameters.txt', 'w') as f:
        for Key in Parameters:
            Value = int_2_str(Parameters[Key])
            f.write(f"{Key}={','.join(Value) if isinstance(Value, list) else Value}\n")

# 📜 Run an R script
def run_Rscript(Rscript):
    subprocess.call(f"R CMD BATCH {Rscript}", shell=True)

# 📊 Count gene-level read counts using annotation
def count_genewise_readcounts(Project_dir, Tools_dir, Result_dir, genomes, Annotation, Gtf_attr_type, bam_file, gene_count_file):
    for ii in range(len(genomes)):
        genome = genomes[ii]

        if genome == "mm9":
            Sub_path_align = "Alignment/minimap/mm9/"
            gtf = "genome/mm9/Annotation/mm9.refGene.gtf" if Annotation == "refGene" else "genome/mm9/Annotation/gencode.vM1.annotation.gtf"
        elif genome == "mm39":
            Sub_path_align = "Alignment/minimap/mm39/"
            gtf = "genome/mm39/Annotation/refGene.gtf"
        elif genome == "hg38":
            Sub_path_align = "Alignment/minimap/hg38/"
            gtf = "genome/hg38/Annotation/hg38.refGene.gtf"

        Result_dir2 = [x.replace("Sub_path/", Sub_path_align) for x in Result_dir]

        for res_path in Result_dir2:
            createFolder(res_path)
            os.chdir(Tools_dir)

            params = {
                "Project_dir": Project_dir,
                "Path": res_path,
                "Annot_file": gtf,
                "Gtf_attr_type": Gtf_attr_type,
                "bam_file_name": bam_file,
                "gene_count_file": gene_count_file
            }
            write_parameters(params)
            run_Rscript("count_genewise_readcounts.R")
            os.chdir(Project_dir)

# 🔢 Get FASTQ read count by dividing total lines by 4
def get_fastq_read_counts(fastq_dirs, fastq_select):
    for fastq_path in fastq_dirs:
        fastq_file = fastq_path + fastq_select
        out_file = fastq_path + "N_fastq_reads.txt"
        command = f"echo $(zcat {fastq_file} | wc -l)/4 | bc > {out_file}"
        print(command)
        subprocess.call(command, shell=True)

# 🌍 Define working directory paths
Tools_dir = os.getcwd()
Project_dir = "/data/groups/anko/mohamed/oxford_nanopore/"

####################################### 📂 Define input/output directories ########################################

fastq_dir = [
    "Data/Nanopore/4_noUV_transcriptome/Mauro/TotalRNA_cDNA_no_UV_rep3/merge_dorado_outputs/",
    "Data/Nanopore/3_SRSF3_DIRCLIP/Mauro/CLIP_cDNA_rep1/merge_dorado_outputs/",
    "Data/Nanopore/3_SRSF3_DIRCLIP/Mauro/CLIP_cDNA_rep2/merge_dorado_outputs/"
]

Result_dir = [
    "Results/Nanopore/4_noUV_transcriptome/Sub_path/Mauro/TotalRNA_cDNA_no_UV_rep3/",
    "Results/Nanopore/3_SRSF3_DIRCLIP/Sub_path/Mauro/CLIP_cDNA_rep1/",
    "Results/Nanopore/3_SRSF3_DIRCLIP/Sub_path/Mauro/CLIP_cDNA_rep2/"
]

# 📌 Set working directory for the project
os.chdir(Project_dir)

# 🔬 Specify genome for alignment and annotation (mm9, mm39, hg38)
genomes = ["mm9"]

################################### 📊 Step 1: Quality Check with FastQC & MultiQC ####################################

Sub_path_fqc = "fastqc/"
fastq_select = "All_merged.fastq.gz"  # You can switch to .fastq if needed

# ⛓️ Update result directory paths for fastqc output
Result_dir1 = [x.replace("Sub_path/", Sub_path_fqc) for x in Result_dir]

# 🚀 Run FastQC across samples
run_fastqc(fastq_dir, fastq_select, Result_dir1, threads=5)

################################### 🔢 Step 2: Count FASTQ Read Depths ###############################################

get_fastq_read_counts(fastq_dir, fastq_select)

################################### 🧬 Step 3: Sequence Alignment with minimap2 ######################################

bam_file = "sorted.bam"
bed_file = "sorted.bed"

run_minimap2_aligner(fastq_dir, fastq_select, genomes, Result_dir, bam_file, bed_file)

################################### 📈 Step 4a: Gene-level Read Count Estimation #####################################

Annotation = "Gencode"  # Alternative: "refGene"
Gtf_attr_type = "gene_name"

gene_count_file = f"{bam_file}_{Annotation}_Genewise_readcounts.tsv"

count_genewise_readcounts(Project_dir, Tools_dir, Result_dir, genomes, Annotation, Gtf_attr_type, bam_file, gene_count_file)

################################### 📈 Step 4b: Gene-level Counts from Filtered BAM ##################################

bam_file = "filt_q8.bam"
gene_count_file = f"{bam_file}_{Annotation}_Genewise_readcounts.tsv"

count_genewise_readcounts(Project_dir, Tools_dir, Result_dir, genomes, Annotation, Gtf_attr_type, bam_file, gene_count_file)

########################################################################################################################








