#!/bin/bash
#SBATCH --job-name=NF_rnafusion
#SBATCH --time=10-00:00:00
#SBATCH --partition=longq
#SBATCH --qos=longq
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=lab_cresswell
#SBATCH --out=/nobackup/lab_cresswell/vgazziero/logs/GBM_fusion.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/GBM_fusion.err


# Apptainer cache directory
NXF_SINGULARITY_CACHEDIR=/nobackup/lab_cresswell/vgazziero/nextflow_pipelines/singularity_images
export NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR

# java virtual machines are started with 1 GB of mem and can reach a maximum of 64 GB
NXF_OPTS='-Xms1g -Xmx64g'
export NXF_OPTS=$NXF_OPTS

# nextflow config
nextflow_config=/nobackup/lab_cresswell/vgazziero/nextflow_pipelines/cemm.config
base=/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/rnaseq
# Load CeMM modules
module load Nextflow/26.04.0|| exit 1
module load apptainer/1.1.9 || exit 1

export NXF_SYNTAX_PARSER=v1

nextflow run nf-core/rnafusion \
  -r 4.1.3 \
  -profile singularity \
  --tools arriba,ctatsplicing,fusioncatcher,starfusion,stringtie,fastp,salmon,fusioninspector,fusionreport \
  --input ${base}/input.csv \
  --cosmic_username virginiaanna.gazziero@phd.units.it --cosmic_passwd fotografia%Orolog1overd3 \
  --genomes_base /nobackup/lab_cresswell/vgazziero/synteny_prj/ewing_sarcoma/fusion_calling/references \
  --outdir ${base}/results \
  -c ${nextflow_config} \
  --skip_vis FALSE
