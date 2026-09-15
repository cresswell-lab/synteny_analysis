library(tidyverse)

rm(list=ls())

path = '/nobackup/lab_cresswell/vgazziero/data/GBM'

metadata = list.files(path, pattern = 'tsv', full.names = T)
metadata = read.table(metadata, sep = '\t', header = T)
chip_input = metadata %>% 
  filter(library_strategy == 'ChIP-Seq') %>% 
  filter(str_detect(sample_title, pattern = 'input|Input|INPUT')) %>% 
  tidyr::separate_rows(fastq_ftp, sep = ';')

x = chip_input %>% 
  pull(run_accession) %>% 
  unique %>% 
  paste(collapse = '|')

fastq_1 = list.files(path, pattern = x, full.names = T) %>% grep('1.fastq.gz', ., value = T)
fastq_2 = list.files(path, pattern = x, full.names = T) %>% grep('2.fastq.gz', ., value = T)

f1 = tibble(
  fastq_1 = fastq_1
) %>% 
  mutate(sample = str_extract(fastq_1, x))
f2 = tibble(
  fastq_2 = fastq_2
) %>% 
  mutate(sample = str_extract(fastq_2, x))

fastqs = full_join(f1, f2)

input = metadata %>% 
  filter(library_strategy == 'ChIP-Seq') %>% 
  filter(str_detect(sample_title, pattern = 'input|Input|INPUT')) %>% 
  full_join(., fastqs, by = join_by('run_accession' == 'sample')) %>% 
  dplyr::select(sample_title, study_accession,sample_alias, fastq_1, fastq_2) %>% 
  mutate(control_type = 'input') %>% 
  dplyr::rename(lab = study_accession)
  
input = input %>% 
  tidyr::separate(sample_title, into = 'cell_line', sep = '-', remove = F) %>% 
  mutate(cell_line = gsub(' Input', '', cell_line)) %>% 
  mutate(cell_line = gsub(' ', '_', cell_line)) %>% 
  relocate(sample_alias, .before = sample_title) %>% 
  relocate(lab, .after = sample_alias) %>% 
  mutate(sample_title = NULL) %>% 
  relocate(control_type, .before = fastq_1)
  
write.csv(input, '/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/input.csv', quote = F, row.names = F)



