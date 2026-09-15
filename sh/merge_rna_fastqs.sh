#!/bin/bash
#SBATCH --job-name=merge_fastq
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shortq 
#SBATCH --qos=shortq 
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=/nobackup/lab_cresswell/vgazziero/logs/merge_fastq.out
#SBATCH --error=/nobackup/lab_cresswell/vgazziero/logs/merge_fastq.err

#!/bin/bash

csv="/nobackup/lab_cresswell/vgazziero/PRJNA674799_A673_RNAseq/merge_fastq.csv"
fastq_dir="/nobackup/lab_cresswell/vgazziero/PRJNA674799_A673_RNAseq"      # adjust: where the raw fastq files live
outdir="/nobackup/lab_cresswell/vgazziero/PRJNA674799_A673_RNAseq/merged_fastqs"        # adjust: where merged files should go
mkdir -p "$outdir"

samples=$(awk -F',' '{print $2}' "$csv" | sort -u)

for s in $samples; do
  # pull SRA accessions for this sample into a bash array
  mapfile -t sra_acc < <(awk -F',' -v samp="$s" '$2==samp {print $1}' "$csv")

  if [ "${#sra_acc[@]}" -ne 2 ]; then
  echo "WARNING: sample $s has ${#sra_acc[@]} accessions, expected 2 — skipping" >&2
  continue
  fi

  echo "Merging sample $s: ${sra_acc[0]} + ${sra_acc[1]}"

  # concatenate R1 files
  cat "${fastq_dir}/${sra_acc[0]}_1.fastq.gz" \
  "${fastq_dir}/${sra_acc[1]}_1.fastq.gz" \
  > "${outdir}/${s}_1.fastq.gz"

  # concatenate R2 files
  cat "${fastq_dir}/${sra_acc[0]}_2.fastq.gz" \
  "${fastq_dir}/${sra_acc[1]}_2.fastq.gz" \
  > "${outdir}/${s}_2.fastq.gz"
done