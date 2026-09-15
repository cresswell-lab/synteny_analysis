#!/usr/bin/env bash
#SBATCH --job-name=hic2cool
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=mediumq 
#SBATCH --qos=mediumq 
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/hic2cool_%a.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/hic2cool_%a.err
#SBATCH --array=0-7

source ~/.bashrc
conda activate hic2cool

path=/nobackup/lab_cresswell/vgazziero/data/GBM/matrices/hic
files=$(find $path -type f -name "*.hic")

mapfile -t files < <(find "$path" -type f -name "*.hic")

mcool_files=/nobackup/lab_cresswell/vgazziero/data/GBM/matrices/mcool
mkdir -p $mcool_files

i=$SLURM_ARRAY_TASK_ID
echo $i

file="${files[$i]}"
echo $file

new_id=${file%.*}.mcool
new_id=$mcool_files/$(basename "$new_id")

echo $new_id

hic2cool convert $file $new_id

