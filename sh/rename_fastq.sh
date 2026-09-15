#!/bin/bash
#SBATCH --job-name=renaming_reads
#SBATCH --array=1-28
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shortq 
#SBATCH --qos=shortq 
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/rename_fastq_%A_%a.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/rename_fastq_%A_%a.err

set -euo pipefail

csv="/nobackup/lab_cresswell/vgazziero/shallow_sarek/fix_fastq.csv"

# pull just this task's line
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$csv")
IFS=',' read -r fastq old_id new_id <<< "$line"

echo "Task $SLURM_ARRAY_TASK_ID: $fastq  ($old_id -> $new_id)"

tmp="${fastq%.fastq.gz}_renamed.fastq.gz"

zcat "$fastq" | sed "s/$old_id/$new_id/g" | gzip > "$tmp"
