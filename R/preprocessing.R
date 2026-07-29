# prepare input sample sheet
rm(list=ls())
library(tidyverse)

path = '/nobackup/lab_cresswell/vgazziero'
cut_and_run = list.files(path, pattern = 'CutRun', full.names = T)

# slicing the metadata to get only the files i am interested, one path at time

metadata = lapply(cut_and_run, function(pp) {
  metadata = list.files(pp, pattern = 'ena_metadata.tsv', full.names = T) %>% 
    read.table(., header = T, sep = '\t') %>% 
    filter(str_detect(pattern = 'IgG', sample_title)) %>%
    separate_rows(sep = ';', fastq_ftp) %>%
    dplyr::select(-c(submitted_ftp, bam_ftp)) %>%
    separate(fastq_ftp,into = c('rm', 'R'), remove = F, sep = '_') %>%
    mutate(rm = NULL,
           sample_id = paste(run_accession, R, sep = '_'))


  fastqs = tibble(fastq_path = list.files(pp, pattern = 'fastq', full.names = T)) %>%
    mutate(sample_id = gsub(paste0(pp, '/'), '', fastq_path))

  metadata = metadata %>%
    full_join(., fastqs)

  write.table(metadata, file = file.path(pp, 'ena_metadata_corrected.tsv'), sep = '\t', quote = F, col.names = T, row.names = F)
  
  return(metadata)
  
})

# read the correct metadata table 

metadata = lapply(cut_and_run, function(pp) {
  list.files(pp, pattern = 'ena_metadata_corrected.tsv', full.names = T) %>% 
    read.table(., header = T, sep = '\t')
}) %>% 
  bind_rows() %>% 
  filter(!is.na(run_accession)) %>% 
  filter(str_detect(pattern = 'TC71|A673', sample_title))

# creating the input for sarek alignment 
# order must be: 
# patient, sex, status, sample, lane, fastq_1, fastq_2

input_csv = metadata %>% 
  dplyr::select(sample_title, fastq_path) %>% 
  tidyr::separate(sample_title, into = 'patient', sep = '-', remove = F) %>% 
  mutate(sex = NA, 
         status = 1, 
         lane = 1) %>% 
  rename(sample = sample_title) %>% 
  relocate(sample, .after = patient) %>% 
  mutate(fastq_path = gsub('.fastq.gz', '_renamed.fastq.gz', fastq_path)) %>% 
  mutate(fastq = case_when(str_detect(pattern = '_1_renamed.fastq.gz', string = fastq_path) ~ 'fastq_1', 
                           str_detect(pattern = '_2_renamed.fastq.gz', string = fastq_path) ~ 'fastq_2', 
                           .default = NA)) %>% 
  pivot_wider(names_from = fastq, values_from = fastq_path) %>% 
  mutate(sample = gsub('-', '_', sample)) %>% 
  mutate(sample = gsub('_CUT&Tag ', '_', sample)) %>% 
  relocate(sample, .after = status) 
write.table(input_csv, '/nobackup/lab_cresswell/vgazziero/shallow_sarek/input.csv', sep =',', quote = F, row.names = F, col.names = T)
  

# csv to change the ids in the fastqs

fix_fastq = metadata %>% 
  dplyr::select(sample_title, fastq_path, run_accession) %>% 
  mutate(sample = gsub('-', '_', sample_title), 
         sample_title = NULL) %>% 
  mutate(sample = gsub('_CUT&Tag ', '_', sample)) 
write.table(fix_fastq, '/nobackup/lab_cresswell/vgazziero/shallow_sarek/fix_fastq.csv', sep =',', quote = F, row.names = F, col.names = F)


