#!/usr/bin/env bash
#SBATCH --job-name=rename_chr
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=longq 
#SBATCH --qos=longq
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/rename_chr_%A.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/rename_chr_%A.err

# source ~/.bashrc
# conda activate neoloop

module load apptainer

file=$1
echo $file
sample=${file%.*}
filt_file=${sample}_filtered.mcool

sample=$(basename "$sample")

output_path=/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/EagleC/cnv/$sample
mkdir -p $output_path

cd /nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/EagleC/cnv

resolutions=(5000 10000 50000)
i=$SLURM_ARRAY_TASK_ID

res="${resolutions[$i]}"

cna_call=${output_path}/${sample}_${res}.CNV_profile.bedGraph

image=/nobackup/lab_cresswell/vgazziero/nextflow_pipelines/singularity_images/neoloop%3A0.4.3.post2--pyhdfd78af_0
cache=/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/EagleC/cnv/.cache
log_file=$output_path/cnv_call.log

apptainer exec --bind /nobackup:/nobackup "$image" \
   python /home/vgazziero/nobackup/synteny_analysis/sh/hic_data/rename_chr.py $filt_file