rm(list=ls())
library(tidyverse)

path = '/nobackup/lab_cresswell/vgazziero/PRJNA674799_A673_RNAseq'

metadata = read.delim(file.path(path, 'metadata.tsv'), sep = '\t', header = T)
filter_samples = paste0('GSM4885', c(656, 661, 657, 662, 658, 663))
filt_metadata = metadata %>% 
  filter(sample_alias %in% filter_samples)

df_rename = filt_metadata %>% 
  dplyr::select(run_accession, sample_alias)
write.table(df_rename, 'PRJNA674799_A673_RNAseq/merge_fastq.csv', sep = ',',col.names = F, row.names = F, quote = F)

sample_ids = filt_metadata %>% 
  dplyr::select(sample_alias, sample_title) %>% 
  mutate(sample = case_when(
    sample_alias == 'GSM4885658' ~ 'iEF_rescue_R1',
    sample_alias == 'GSM4885661' ~ 'iLuc_empty_R2', 
    sample_alias == 'GSM4885663' ~ 'iEF_rescue_R2',
    sample_alias == 'GSM4885656' ~ 'iLuc_empty_R1',
    sample_alias == 'GSM4885657' ~ 'iEF_empty_R1',
    sample_alias == 'GSM4885662' ~ 'iEF_empty_R2'
  )) %>% 
  dplyr::select(sample_alias, sample) %>% 
  distinct()
  
run_acc_id = filt_metadata %>% 
  pull(sample_alias) %>% 
  unique %>% 
  paste(collapse = '|')

fastqs_1 = list.files(file.path(path, 'merged_fastqs'), pattern = run_acc_id, full.names = T) %>% grep('1.fastq.gz', ., value = T)
fastqs_2 = list.files(file.path(path, 'merged_fastqs'), pattern = run_acc_id, full.names = T) %>% grep('2.fastq.gz', ., value = T)

id = gsub(pattern = '/nobackup/lab_cresswell/vgazziero/PRJNA674799_A673_RNAseq/merged_fastqs/', '',fastqs_1) %>% 
  gsub('_1.fastq.gz', '', .)
input = tibble(
  id = id,
  fastq_1 = fastqs_1, 
  fastq_2 = fastqs_2, 
  strandedness = 'forward'
) %>% 
  full_join(., sample_ids, by = join_by('id' == 'sample_alias')) %>% 
  dplyr::select(-id) %>% 
  relocate(sample, .before = fastq_1)
write.table(input, 'synteny_prj/fusion_calling/input.csv', sep = ',', quote = F, row.names = F, col.names = T)


path = '/nobackup/lab_cresswell/vgazziero/data/GBM'
res = 'synteny_prj/GBM/rnaseq'
dir.create(res)

metadata = read.delim(file.path(path, 'metadata.tsv'), sep = '\t', header = T)
filt_metadata = metadata %>% 
  filter(library_strategy == 'RNA-Seq')

sra = filt_metadata$run_accession %>% paste(collapse = '|')
fastq = list.files(path, pattern = sra, full.names = T)
setdiff(input$fastq_1, fastq)


input = filt_metadata %>% 
  dplyr::select(run_accession, sample_title) %>% 
  mutate(fastq_1 = paste0(path, '/', run_accession, '.fastq.gz'), 
         run_accession = NULL, 
         fastq_2 = '', 
         strandedness = 'forward') %>% 
  mutate(sample_title = gsub('-', '_', sample_title)) %>% 
  mutate(sample_title = gsub('RNAseq_', '', sample_title)) %>% 
  rename(sample = sample_title)
input

write.csv(input, file.path(res, 'input.csv'), quote = F, row.names = F, col.names = T)