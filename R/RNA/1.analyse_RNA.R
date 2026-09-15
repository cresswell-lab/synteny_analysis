rm(list=ls())
setwd('/research/lab_cresswell/vgazziero/synteny_prj')
library(tidyverse)
library('DESeq2')
library(org.Hs.eg.db)
library(clusterProfiler)
source('synteny_prj/scripts/utils.R')

res_path = '/research/lab_cresswell/vgazziero/synteny_prj/r_obj'
omic = 'RNAseq'
res_path = file.path(res_path, omic)
dir.create(res_path, recursive = T)

img_path = '/research/lab_cresswell/vgazziero/synteny_prj/img'
img_path = file.path(img_path, omic)
dir.create(img_path, recursive = T)

# accounting for "normal" RNA counts. then try to download the fastqs and run STAR to detect fusions

counts_path = '/research/lab_cresswell/vgazziero/synteny_prj/datasets/ewing_sarcoma/rnaseq_counts/rnaseq_counts.tsv'
rna_counts = read.table(counts_path, sep = '\t', stringsAsFactors = FALSE)
colnames(rna_counts) = rna_counts[1,]
rna_counts= rna_counts[-1,]

# map entrez to gene symbols
entrez_ids = rna_counts[,'GeneID']

genes = convert_names(entrez_ids, from = 'ENTREZID', to = 'SYMBOL', org = org.Hs.eg.db)

rna_geo = getGEO('GSE185130')
metadata = rna_geo$GSE185130_series_matrix.txt.gz@phenoData@data %>%
  dplyr::select(geo_accession, characteristics_ch1) %>%
  mutate(condition = ifelse(str_detect(characteristics_ch1, 'iLuc'), 'EFEndo', 'EFKD')) %>% 
  dplyr::select(condition) %>% 
  mutate(condition = factor(condition, levels = c('EFEndo', 'EFKD')))

# remove duplicated entrez
genes = rna_counts %>% 
  as_tibble() %>% 
  left_join(., genes, by = join_by('GeneID' == 'ENTREZID')) %>% 
  dplyr::select(GeneID, SYMBOL) %>% 
  distinct()

rna_counts = rna_counts %>% 
  as_tibble() %>% 
  # left_join(., genes, by = join_by('GeneID' == 'ENTREZID')) %>%
  mutate_all(., function(x) as.numeric(as.character(x))) %>% 
  # relocate(SYMBOL, .before = GeneID) %>% 
  # dplyr::select(-GeneID) %>% 
  tibble::column_to_rownames('GeneID') %>% 
  as.matrix()

# create deseq2 obj and perform dge 
dds = launch_deseq(raw_counts = rna_counts, 
                   coldata = metadata, 
                   design = '~ condition', 
                   genes = genes$SYMBOL, 
                   ref_lv = 'EFEndo', 
                   sGS = 3, 
                   minCounts = 10)
saveRDS(dds, file = file.path(res_path, 'dds.rds'))

# run DE
res <- results(dds, contrast=c("condition","EFKD","EFEndo"), tidy = T)
res = left_join(res, genes, by = join_by(row == GeneID))
# write.table(res, 'cu_vs_atn/analysis/res/cuR_vs_cu_results_v2.csv', sep = ',', quote = F, row.names = F, col.names = T)
saveRDS(res, file.path(res_path, 'de_res.rds'))

# pca check (note that four samples are really few...)
vsd_res = vst(dds, blind=FALSE)
plotPCA(vsd_res, intgroup="condition") +
  theme_bw() +
  # scale_color_manual(values = c(Cu = "#009FBD", ATN = "#EE4266"))+ 
  ggrepel::geom_text_repel(aes(label = name)) +
  # geom_text(aes(label = name)) +
  geom_point()

# normalise data using vst
res_norm = assay(vst(dds, blind = FALSE))
res_norm = res_norm %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('gene_id') %>% 
  left_join(.,  genes, by = join_by(gene_id == GeneID)) %>% 
  rename(gene_name = SYMBOL) %>% 
  dplyr::relocate(gene_name, .after = gene_id)
saveRDS(res_norm, file.path(res_path, 'vst_expression.rds'))

