library(ACE)
library(tidyverse)
library(QDNAseq.hg38)
library(copynumber)

# remotes::install_github("aroneklund/copynumber")
# devtools::install_github("asntech/QDNAseq.hg38@main")
# BiocManager::install('ACE')

bam_path = '/nobackup/lab_cresswell/vgazziero/cut_and_tag_processing/cna_calling'
bams = list.files(bam_path, full.names = T, pattern = 'bam$')

out_dir = file.path(bam_path, 'ACE_results_new_alignment')
dir.create(out_dir, showWarnings = F)

bin_size = 1000

runACE(bam_path, filetype ='bam', binsizes = bin_size, ploidies = c(2), imagetype='png', genome = "hg38", outputdir = out_dir)

object = readRDS(paste0(out_dir,"/",bin_size,"kbp.rds"))
object@phenoData@data$name
object@phenoData@data

# Check number of reads
number_of_reads <- data.frame(Sample = object@phenoData@data$name, TotalReads = object@phenoData@data$total.reads, UsedReads = object@phenoData@data$used.reads)

# Subset for samples of interest

A673_samples <- grep("A673" , colnames(object)) #here you can subset for whatever samples you want, e.g. grep("cell_line_name")
#samples_subset <- 1:72 #if you want to look at all samples

# Get the bins

results_list = lapply(A673_samples, function(i) {
    sample_id_index = i
    sol_number      = 1 #taking the "best" solution
  
    sq_model = squaremodel(object, QDNAseqobjectsample = sample_id_index, 
                         cellularities = seq(95, 100, by = 0.1), 
                         ptop = 2.5, pbottom = 1.5, prows = 20)
  
    temp_results = objectsampletotemplate(object, index = sample_id_index)
  
  # Add column to identify the sample
  temp_results$Sample <- colnames(object)[sample_id_index]
  return(temp_results)
})

names(results_list) = grep("A673" , colnames(object), value = T) 

# Set parameters
sample_id <- "A673" #cell line
cna_mat = NULL
psit_mat = NULL
rho_mat = NULL
model_check = NULL
sol_number = 1


for (sample in names(results_list)) {
  single_purity <- results_list[[sample]] #extract dataframe for each purity
  
  # Calculate logged copynumbers
  single_purity$log2_copynumber <- log(abs(single_purity$copynumber), base = 2)
  
  # Add a position column
  single_purity$pos <- (single_purity$end + single_purity$start)/2
  
  # Get rid of all the columns we do not use
  single_purity <- single_purity %>%
  dplyr::select(chr, pos, Sample, log2_copynumber) %>%
  filter(chr != "Y")
  
  # Group to get the log2value of each sample at each position of each chromosome
  single_purity <- single_purity %>%
  group_by(chr, pos, Sample) %>%
  dplyr::slice(1) %>%  
  ungroup()
  
  single_purity$chr <- as.numeric(single_purity$chr)
  single_purity <- single_purity %>% arrange(chr, pos)
  
  # change format so that each sample is one column
  single_purity_wide <- single_purity %>%
  pivot_wider(
    names_from = Sample,
    values_from = log2_copynumber
  )
  
  # exclude NAs
  single_purity_wide <- na.omit(single_purity_wide)
  
  single_purity_wide <- as.data.frame(single_purity_wide)
  head(single_purity_wide)
  
  res_segs = pcf(single_purity_wide, assembly = "hg38", gamma = 40)
  
  # Data preparation
  # res_segs contains the mean logged copynumbers per segment --> need the data not logged
  res_segs[,7] <- 2^(res_segs[,7])
  
  # single_purity_wide contains the log2 ratio for each sample per bin --> also need it not logged for the function
  res_bins <- single_purity_wide
  res_bins[, c(3:ncol(res_bins))] <- 2^(single_purity_wide[, c(3:ncol(res_bins))])
  
  # Remove the first column and make the sample name the colname (instead of mean)
  sample_label <- res_segs[[1]]
  colnames(res_segs)[7] <- sample_label[1]
  res_segs <- res_segs[,-1]
  
  # Plot
  # Build template dataframe
  template <- data.frame(
    bin = c(1:nrow(res_bins)),
    chr = rep(res_segs$chrom, times = res_segs$n.probes),
    start = rep(res_segs$start.pos, times = res_segs$n.probes),
    end = rep(res_segs$end.pos, times = res_segs$n.probes),
    copynumbers = res_bins[[sample]],
    segments = rep(res_segs[[sample]], times = res_segs$n.probes)
  )
    
  # Run squaremodel
  result_sqm <- squaremodel(template, cellularities = seq(95, 100, by = 0.1), ptop = 2.5, pbottom = 1.5, prows = 20)
  #ploidy 1.5-2.5 for diploid, 2.5-3.5 for triploid and 3.5-4.5 for tetraploid
    
  # Include number of reads
  used_reads <- number_of_reads$UsedReads[number_of_reads$Sample == sample]
    
  # Get adjusted segments
  segmentdf = getadjustedsegments(template, cellularity = result_sqm$minimadf[sol_number,"cellularity"], 
                                  ploidy = result_sqm$minimadf[sol_number,"ploidy"])
  psit  = result_sqm$minimadf[sol_number,"ploidy"]
  rho   = result_sqm$minimadf[sol_number,"cellularity"]
  segmentdf$Chromosome <- as.numeric(segmentdf$Chromosome)
  segmentdf <- segmentdf %>% arrange(Chromosome, Start)
    
  # Plot
  #Color scheme
  cols = c(c("0" = "#1981be", "1" = "#56B4E9", "2" = "grey", "3" = "#E69F00", "4" = "#ffc342",
             "5" = "#FFAA42", "6" = "#FF9142", "7" = "#FF7742", "8" = "#FF5E42"), 
           rep("#FF4542", times = 1000 - 8))
  names(cols)[(9:1000)+1] = 9:1000
    
  # Prepare df for plotting
  plt.df = data.frame(genome.bin = 1:nrow(na.omit(template)), 
                    chr = factor(paste0("chr",rep(segmentdf$Chromosome,
                                  times = segmentdf$Num_Bins)),
                                  levels = paste0("chr",c(1:22))), 
                    start = na.omit(template)$start,
                    end = na.omit(template)$end,
                    Log2ratio = log(na.omit(template)$copynumbers, base = 2),
                    Call = as.character(rep(segmentdf$Copies, times = segmentdf$Num_Bins)),
                    Segment = log(na.omit(template)$segments, base = 2))
  
  plt.df$Call[plt.df$Call == "-1"] <- NA
  cn_levels <- unique(plt.df$Call)
  cn_levels <- cn_levels[order(as.numeric(cn_levels))]
  plt.df$Call <- factor(plt.df$Call, levels = cn_levels)
  cn_levels <- cn_levels[cn_levels != "NA"]

    
  # Make the plot
  plot = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
  geom_hline(yintercept = c(-3,-2,-1,0,1,2,3), lty = c("solid"), lwd = 0.2) +
  geom_point() +
  ylab("log2ratio") +
  scale_colour_manual(limits = cn_levels, values = cols[cn_levels]) +
  scale_x_continuous(name = "Chromosomes", labels = c(1:22), 
                      breaks = as.vector(c(1, cumsum(table(plt.df$chr))[-22]) + 
                                            (table(plt.df$chr) / 2))) + 
  geom_vline(xintercept = c(1, cumsum(table(plt.df$chr))), lty = "dotted") +
  ggtitle(sample_id) + 
  scale_y_continuous(limits=c(-2,2), oob=scales::squish) + 
  geom_point(aes(y = Segment), color="#000000", size = 0.1) +
  #theme_cowplot() +
  labs(subtitle = paste0("ploidy=",
                   result_sqm$minimadf[sol_number,"ploidy"], ", purity=",
                   result_sqm$minimadf[sol_number,"cellularity"],", used reads=",used_reads)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 26),
        plot.subtitle = element_text(hjust = 0.5, size = 22),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
  
  print(plot)

  dir.create(file.path(out_dir, "img"))
  ggsave(file.path(out_dir, "img", paste0(sample, '_', bin_size, '_kbp_ACE_cn_profile.png')), 
        plot = plot, height = 4, width = 15, dpi = 300, units = 'in')
  #ggsave(paste0(out_dir,"/plots_thesis/",sample,"_",
              #bin_size,"kbp_ACE_cn_profile.png"), 
        #plot = plot,
        #height = 3.5, width = 15)
    
  # Collect calls
  plt.df$Call <- as.character(rep(segmentdf$Copies, times = segmentdf$Num_Bins)) 
  
  cna_mat = cbind(cna_mat, plt.df$Call) 
  psit_mat = cbind(psit_mat, psit)
  rho_mat = cbind(rho_mat, rho)
  model_check = cbind(model_check, result_sqm)
  
  # saveRDS(file = paste0(out_dir, "/RDS_thesis/STA-NB-13/", sample, "_segments.rds"), psit_mat)
}

# Add sample names to the calls
colnames(cna_mat) <- sample_names
colnames(psit_mat) <- sample_names
colnames(rho_mat) <- sample_names



