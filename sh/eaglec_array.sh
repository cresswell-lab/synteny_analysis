#!/usr/bin/env bash

path=/nobackup/lab_cresswell/vgazziero/data/GBM/matrices/mcool
files=$(find $path -type f -name "*.mcool")

mapfile -t files < <(find "$path" -type f -name "*.mcool")

for f in "${files[@]}"; do
    
    echo $f
    sbatch /home/vgazziero/nobackup/synteny_analysis/sh/run_eagleC_raw_matrices.sh $f

done
