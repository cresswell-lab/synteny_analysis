# utils 
convert_names = function(x, from, to, org) {
  
  # convert
  gene_df <- bitr(
    x,
    fromType = from,
    toType   = to,
    OrgDb    = org
  )  
  
  # x <- gene_df[,to]
  
  return(gene_df)
}

# deseq2 utils

read_raw_counts = function(path) {
  
  data = read.table(path, sep = "\t", header = T)
  
  gene_names = data %>% 
    dplyr::select(gene_id, gene_name)
  
  raw_counts = data %>% 
    dplyr::select(-gene_name) %>% 
    tibble::column_to_rownames('gene_id') %>% 
    mutate_all(~ as.integer(.)) %>% 
    as.matrix()
  
  tt = list('raw_counts' = raw_counts, 
            'map_genes_ensembl' = gene_names) 
  return(tt)
}

launch_deseq = function(raw_counts, 
                        coldata, 
                        design, 
                        genes, 
                        ref_lv, 
                        sGS = 3, 
                        minCounts = 10) {
  
  # construct the deseqdataset object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                        colData = coldata,
                                        design = as.formula(design))
  
  # add metadata features to the deseq object
  mcols(dds) <- DataFrame(mcols(dds), genes)
  # mcols(dds)
  
  # specify the factor levels
  dds$condition <- relevel(dds$condition, ref = ref_lv)
  
  smallestGroupSize <- sGS
  keep <- rowSums(counts(dds) >= minCounts) >= smallestGroupSize
  dds <- dds[keep,]
  
  # launch the differential expression 
  dds_v2 <- DESeq(dds)
  
  return(dds_v2)
}


