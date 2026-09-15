rm(list=ls())
library(ACE)
library(tidyverse)
library(QDNAseq.hg38)
library(copynumber)
library(patchwork)
source('/nobackup/lab_cresswell/vgazziero/synteny_analysis/R/colors.R')
source('/nobackup/lab_cresswell/vgazziero/synteny_analysis/R/CNA/utils.R')
# source('/nobackup/lab_cresswell/vgazziero/synteny_analysis/R/utils/runACE.R')

if(run_ace) {
  bam_path = '/nobackup/lab_cresswell/vgazziero/cut_and_tag_processing/cna_calling'
  bams = list.files(bam_path, full.names = T, pattern = 'bam$')
  
  out_dir = file.path(bam_path, 'ACE_results_new_alignment')
  dir.create(out_dir, showWarnings = F)
  
  runACE(bam_path, filetype ='bam', binsizes = bin_size, ploidies = c(2), imagetype='png', genome = "hg38", outputdir = out_dir)
}

# setting paths 
res_path = '/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/cna'
dir.create(res_path, showWarnings = F)
rds_path = file.path(res_path, 'rds')
dir.create(rds_path, showWarnings = F)
img_path = file.path(res_path, 'img')
dir.create(img_path, showWarnings = F)
path = '/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/chipseq/results/input_samples'
cell_lines = read.table('/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/chipseq/input.csv', sep = ',', header = T) %>% 
  pull(cell_line) %>% 
  unique

# load data
ace_res = lapply(cell_lines, function(cl) {
  print(cl)
  pp = file.path(path, cl, 'cnv', 'ace')
  pp = list.files(pp, pattern = 'ace_out_', full.names = T)
  df = file.path(pp, paste0(bin_size, 'kbp.rds'))
  readRDS(df)
})
names(ace_res) = cell_lines

# get number of reads
number_of_reads = lapply(ace_res, function(x) {
  data.frame(Sample = x@phenoData@data$name, TotalReads = x@phenoData@data$total.reads, UsedReads = x@phenoData@data$used.reads)
}) %>% 
  bind_rows()

# Get the bins
results_list = lapply(ace_res, function(x) {
  
  id = 1
  
  sq_model = squaremodel(x, QDNAseqobjectsample = id, 
                         cellularities = seq(95, 100, by = 0.1), 
                         ptop = 3.5, pbottom = 1.5, prows = 20)
  
  temp_results = objectsampletotemplate(x, index = id)
  
  # Add column to identify the sample
  temp_results$Sample <- colnames(x)[id]
  return(temp_results)
}) 


res = lapply(names(results_list), function(id) {
  
  ploidy_borders = exp_psi %>% 
    filter(cell_line == id) %>% 
    mutate(top = ploidy+.5, 
           bottom = ploidy-.5)
  
  x = results_list[[id]]
  print(id)
  ploidy_borders
  
  extract_results(x, ptop = ploidy_borders$top, pbottom = ploidy_borders$bottom, number_of_reads = number_of_reads, cols = cols, id = id)
})
names(res) = names(results_list)
saveRDS(res, file.path(rds_path, 'squaremodel_results.rds'))

plts = lapply(res, function(x) {x$plot})
patchwork::wrap_plots(plts)

lapply(names(plts), function(x) {
  ggsave(plot = plts[[x]], filename = file.path(img_path, paste0(x, '_cna_profile.png')), width = 15, height = 4, dpi = 300, units = 'in')
  ggsave(plot = plts[[x]], filename = file.path(img_path, paste0(x, '_cna_profile.pdf')), width = 15, height = 4, dpi = 300, units = 'in')
})

