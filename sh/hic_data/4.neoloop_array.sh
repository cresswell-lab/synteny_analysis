#!/usr/bin/env bash

path=/nobackup/lab_cresswell/vgazziero/data/GBM/matrices/mcool
files=$(find $path -type f -name "*_filtered.mcool")

mapfile -t files < <(find "$path" -type f -name "*_filtered.mcool")

for f in "${files[@]}"; do
    
    echo $f
    sbatch /home/vgazziero/nobackup/synteny_analysis/sh/hic_data/array_calculate_cnv_hic.sh $f

done
