#!/bin/bash
#SBATCH --job-name=run_ACE
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shortq 
#SBATCH --qos=shortq 
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/run_ACE_%j.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/run_ACE_%j.err

conda run -n ACE Rscript /nobackup/lab_cresswell/vgazziero/synteny_analysis/R/CNA/1.1.joint_segmentation.R
