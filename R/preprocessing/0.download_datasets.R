# install.packages('tidyverse')
# install.packages("BiocManager")
# BiocManager::install(version = "3.22")
# BiocManager::install('ComplexHeatmap')
# install.packages('remotes')
# BiocManager::install('GEOquery')
library(tidyverse)
library(GEOquery)

# try to fetch data of Erwing sarcoma 

# datasets = c('GSE185125', 
#              'GSE185126', 
#              'GSE185127', 
#              'GSE185128',
#              'GSE185130', 
#              'GSE185131')

dataset = 'GSE162976'
getGEOSuppFiles(GEO = dataset, 
                baseDir = '/nobackup/lab_cresswell/vgazziero/data/GBM/matrices', 
                makeDirectory = F, fetch_files = T)

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162976/suppl//GSE162976_RAW.tar', 
              destfile = '/nobackup/lab_cresswell/vgazziero/data/GBM/matrices/raw_data.tar', 
              method = 'libcurl', 
              mode = 'wb')

# hic, rnaseq, cut and tag and 4C data
lapply(datasets, function(dd) {
  getGEOSuppFiles(dd, baseDir = 'synteny_prj/datasets/ewing_sarcoma', makeDirectory = T)
})

# rnaseq counts
download_rna = FALSE
if(download_rna){
  download.file('https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE185130&format=file&file=GSE185130_raw_counts_GRCh38.p13_NCBI.tsv.gz', 
                destfile = 'synteny_prj/datasets/ewing_sarcoma/rnaseq_counts.tsv.gz', 
                method = "libcurl",
                mode = "wb"
  )
  unzip('synteny_prj/datasets/ewing_sarcoma/rnaseq_counts.tsv.gz')
}



