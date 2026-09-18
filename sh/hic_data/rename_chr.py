#!/usr/bin/env python

import cooler
import sys
import os
import glob

path='/nobackup/lab_cresswell/vgazziero/data/GBM/matrices/mcool'
os.chdir(path)
files=glob.glob('*_filtered.mcool')

for in_cool in files: 
    print(in_cool)
    resolutions = cooler.fileops.list_coolers(in_cool)  # e.g. ['/resolutions/5000', ...]
    for res in resolutions:
        uri = f'{in_cool}::{res}'
        clr = cooler.Cooler(uri)
        chromlist = clr.chromnames
        # print(chromlist)
        chrom_map = {c:'chr'+c for c in chromlist}
        cooler.rename_chroms(clr, chrom_map)

