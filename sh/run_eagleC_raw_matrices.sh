#!/usr/bin/env bash
#SBATCH --job-name=EagleC_parallel_raw
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=mediumq 
#SBATCH --qos=mediumq 
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/EagleC_%J_%a.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/EagleC_%J_%a.err
#SBATCH --array=0-15

source ~/.bashrc
conda activate eagleC

file=$1
echo $file
sample=${file%.*}
sample=$(basename "$sample")

output_path=/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/EagleC/raw/$sample
mkdir -p $output_path

cd $output_path

predictSV --hic-5k $file::/resolutions/5000 \
    --hic-10k $file::/resolutions/10000 \
    --hic-50k $file::/resolutions/50000 \
    -O $sample -g hg38 --balance-type Raw \
    --output-format full