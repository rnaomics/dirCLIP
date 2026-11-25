#!/bin/bash
# submit_all.sh: launch the dirCLIP FAST5 signal pipeline on the HPC cluster.
#   stage 1 (array 0-8): extract FAST5 -> nanopolish index -> eventalign -> NanopolishComp collapse
#   stage 2 (array 0-1): nanocompore sampcomp  (runs only afterok of ALL stage-1 tasks)
# Run from: $BASE/analysis/eventalign
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

EA_JID=$(sbatch --parsable 00_eventalign_array.slurm)
echo "eventalign array submitted: job $EA_JID  (tasks 0-8)"

NC_JID=$(sbatch --parsable --dependency=afterok:${EA_JID} 10_nanocompore_array.slurm)
echo "nanocompore array submitted: job $NC_JID  (afterok:${EA_JID}, tasks 0-1)"

echo
echo "Monitor:  squeue -u \$USER"
echo "Logs:     $PWD/logs/"
