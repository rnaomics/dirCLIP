# conda activate base
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 2 09:17:41 2022

@author: mohamedb
"""

import os
import subprocess
import pandas as pd
import pdb
import sys

# Uncomment if using external module
# sys.path.append('/data/groups/anko/mohamed/oxford_nanopore/Tools/Nanopore_functions')
# import Nanopore_pipeline as nano

# --------------------------------------------------
# 🗂️ Function to create a folder if it doesn't exist
# --------------------------------------------------
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

# --------------------------------------------------
# 🧬 Function to extract FASTA sequences from BED file
# --------------------------------------------------
def get_fasta_from_bed(genome_file, bed_file, fasta_file):
    cmd = f"bedtools getfasta -fi {genome_file} -bed {bed_file} -s -nameOnly > {fasta_file}"
    print(cmd)
    subprocess.call(cmd, shell=True)

    # Compress to .gz (overwrite if exists)
    cmd = f"gzip -f {fasta_file}"
    print(cmd)
    subprocess.call(cmd, shell=True)

# --------------------------------------------------
# 📁 Define global paths and settings
# --------------------------------------------------
parent_dir = "/data/groups/anko/mohamed/oxford_nanopore/"

Mut_paths_all = [
    "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/"
]

Motif_paths_all = [
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/CLIP_cDNA_merge/"
]

ref_subfolder = [
    "ref_P19_WT_TotalRNA_cDNA_no_UV",
    "ref_TotalRNA_cDNA_UV",
    "ref_P19_WT_TotalRNA_cDNA_no_UV",
    "ref_TotalRNA_cDNA_UV"
]

subfolders1 = ["q_8_Q_8"]
subfolders2 = ["2_Inc_del", "1_No_del"]

genome = "mm9"

# --------------------------------------------------
# 🔧 Threshold settings
# --------------------------------------------------
Coverage_thresholds = [50, 100, 200, 300, 400, 500]
Diff_thresholds = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]

filt_by = "FDR"                # Alternative: "pvalue"
sub_folder_suffix = ""         # Optional: "comp_"
wins = [100]                   # Window sizes for motif search
offset = 0                     # Mutation site offset

# --------------------------------------------------
# 🚀 Generate FASTA files from mutation site BEDs
# --------------------------------------------------
for subfolder1 in subfolders1:
    for subfolder2 in subfolders2:
        for i in range(len(Mut_paths_all)):
            Mut_path = Mut_paths_all[i] + f"{subfolder1}/{subfolder2}/{ref_subfolder[i]}/"
            Motif_path = Motif_paths_all[i] + f"{subfolder1}/{subfolder2}/{ref_subfolder[i]}/"

            # Select genome
            if genome == "mm9":
                genome_file = parent_dir + "genome/mm9/mm9.fa"
            elif genome == "mm39":
                genome_file = parent_dir + "genome/mm39/mm39.fa"

            for j in range(len(Coverage_thresholds)):
                for k in range(len(Diff_thresholds)):
                    for win in wins:
                        if offset == 0:
                            bed_path = parent_dir + f"{Mut_path}{sub_folder_suffix}cov_thresh_{Coverage_thresholds[j]}/{filt_by}/"
                            res_path1 = parent_dir + f"{Motif_path}{sub_folder_suffix}cov_thresh_{Coverage_thresholds[j]}/win_{win}/sites/"
                            res_path2 = parent_dir + f"{Motif_path}{sub_folder_suffix}cov_thresh_{Coverage_thresholds[j]}/win_{win}/regions/"
                            res_path3 = parent_dir + f"{Motif_path}{sub_folder_suffix}cov_thresh_{Coverage_thresholds[j]}/win_{win}/regions_centre/"

                            pos_mut_bed_file = bed_path + f"Pos_mut_thresh_{Diff_thresholds[k]}_w_{win}.bed.gz"
                            mut_regions_bed_file = bed_path + f"Pos_mut_thresh_{Diff_thresholds[k]}_w_{win}_region.bed.gz"
                            mut_regions_centre_bed_file = bed_path + f"Pos_mut_thresh_{Diff_thresholds[k]}_w_{win}_region_centre.bed.gz"

                            sites_fasta_file = res_path1 + f"Pos_mut_thresh_{Diff_thresholds[k]}.fasta"
                            regions_fasta_file = res_path2 + f"Pos_mut_thresh_{Diff_thresholds[k]}.fasta"
                            regions_centre_fasta_file = res_path3 + f"Pos_mut_thresh_{Diff_thresholds[k]}.fasta"

                        elif offset > 0:
                            bed_path = parent_dir + f"{Mut_path}{sub_folder_suffix}cov_thresh_{Coverage_thresholds[j]}/"
                            res_path = parent_dir + f"{Motif_path}{sub_folder_suffix}cov_thresh_{Coverage_thresholds[j]}/offset_{offset}_win_{win}/"
                            pos_mut_bed_file = bed_path + f"Pos_mut_thresh_{Diff_thresholds[k]}_offset_{offset}_w_{win}.bed.gz"
                            pos_fasta_file = res_path + f"Pos_mut_thresh_{Diff_thresholds[k]}.fasta"
                            # neg_mut_bed_file = ...
                            # neg_fasta_file = ...

                        # Create output folders
                        createFolder(res_path1)
                        createFolder(res_path2)
                        createFolder(res_path3)

                        # Extract FASTA sequences
                        get_fasta_from_bed(genome_file, pos_mut_bed_file, sites_fasta_file)
                        get_fasta_from_bed(genome_file, mut_regions_bed_file, regions_fasta_file)
                        get_fasta_from_bed(genome_file, mut_regions_centre_bed_file, regions_centre_fasta_file)
                        # get_fasta_from_bed(genome_file, neg_mut_bed_file, neg_fasta_file)
