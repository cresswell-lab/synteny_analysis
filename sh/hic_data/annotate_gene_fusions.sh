#!/usr/bin/env bash
#SBATCH --job-name=annotate_fusions
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=mediumq 
#SBATCH --qos=mediumq 
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/EagleC_annotate_fusions_%J_%a.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/EagleC_annotate_fusions_%J_%a.err
#SBATCH --array=0-7

source ~/.bashrc
conda activate eagleC

export PYENSEMBL_CACHE_DIR=/nobackup/lab_cresswell/vgazziero/.cache

path=/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/EagleC/raw
cd $path
files=(*)

i=$SLURM_ARRAY_TASK_ID

file="${files[$i]}"

sv_call=${path}/${file}/${file}.CNN_SVs.5K_combined.txt
output=${path}/${file}/${file}_gene_fusions.txt

annotate-gene-fusion --sv-file $sv_call \
                       --output-file $output \
                       --buff-size 10000 --skip-rows 1 --ensembl-release 93 --species human
