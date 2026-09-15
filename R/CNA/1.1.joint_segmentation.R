rm(list=ls())
library(ACE)
library(tidyverse)
library(QDNAseq.hg38)
library(copynumber)
library(gridExtra)
source('/nobackup/lab_cresswell/vgazziero/synteny_analysis/R/colors.R')
source('/nobackup/lab_cresswell/vgazziero/synteny_analysis/R/utils/runACE.R')

res_path = '/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/cna'
dir.create(res_path, showWarnings = F)
img_path = file.path(res_path, 'img')
dir.create(img_path, showWarnings = F)
path = '/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/chipseq/results/input_samples'
cell_lines = read.table('/nobackup/lab_cresswell/vgazziero/synteny_prj/GBM/chipseq/input.csv', sep = ',', header = T) %>% 
  pull(cell_line) %>% 
  unique

ace_res = lapply(cell_lines, function(cl) {
  print(cl)
  pp = file.path(path, cl, 'cnv', 'ace')
  pp = list.files(pp, pattern = 'ace_out_', full.names = T)
  df = file.path(pp, '1000kbp.rds')
  readRDS(df)
})
names(ace_res) = cell_lines

# Set parameters
cna_mat_joint = NULL
psit_mat = NULL
rho_mat = NULL
sol_number = 1

# get number of reads
number_of_reads = lapply(ace_res, function(x) {
  data.frame(Sample = x@phenoData@data$name, TotalReads = x@phenoData@data$total.reads, UsedReads = x@phenoData@data$used.reads)
})

# Joint segmentation to show how copy number calling improves
# Get the  bins
combined_results = lapply(ace_res, function(x) {
  
  id = 1
  
  sq_model = squaremodel(x, QDNAseqobjectsample = id, 
                         cellularities = seq(95, 100, by = 0.1), 
                         ptop = 3, pbottom = 1.5, prows = 20)
  
  temp_results = objectsampletotemplate(x, index = id)
  
  # Add column to identify the sample
  temp_results$Sample <- colnames(x)[id]
  return(temp_results)
  # Combine in dataframe
  # combined_results <- rbind(combined_results, temp_results)
}) 

combined_results = lapply(combined_results, function(x) {
  # Calculate logged copynumbers
  x$log2_copynumber <- log(x$copynumber, base = 2)
  
  # add a position column --> position is normal the mean: (end+start)/2 
  x$pos <- (x$end + x$start)/2
  
  return(x)
})

comb_results_filt = lapply(combined_results, function(x) {
  x %>% 
    dplyr::select(chr, pos, Sample, log2_copynumber) %>%
    filter(chr != "Y") %>% 
    group_by(chr, pos, Sample) %>%
    dplyr::slice(1) %>%  
    ungroup() %>% 
    pivot_wider(
      names_from = Sample,
      values_from = log2_copynumber
    ) %>% 
    na.omit() %>% 
    as.data.frame()
})

# perform multisample segmentation
comb_res_segs = lapply(comb_results_filt, function(x) {
  multipcf(x, assembly = "hg38", gamma = 40)
})
  

# Data preparation
# comb_res_segs contains the mean logged copynumbers per segment --> need the data not logged
comb_res_segs[, c(6:ncol(comb_res_segs))] <- 2^(comb_res_segs[, c(6:ncol(comb_res_segs))])

# need to arrange it by increasing chromosome number
comb_res_segs$chrom <- as.numeric(comb_res_segs$chrom)
comb_res_segs <- comb_res_segs %>% arrange(chrom, start.pos)

comb_res_bins <- comb_results_filt_wide
comb_res_bins[, c(3:ncol(comb_res_bins))] <- 2^(comb_results_filt_wide[, c(3:ncol(comb_results_filt_wide))])
comb_res_bins$chr <- as.numeric(comb_res_bins$chr)
comb_res_bins <- comb_res_bins %>% arrange(chr, pos)

# loop over each sample from the patient
sample_names <- colnames(comb_res_segs)[6:ncol(comb_res_segs)]

joint_seg_res = lapply(sample_names, function(s) {
  # Build template dataframe
  template_joint <- data.frame(
    bin = c(1:nrow(comb_res_bins)),
    chr = rep(comb_res_segs$chrom, times = comb_res_segs$n.probes),
    start = rep(comb_res_segs$start.pos, times = comb_res_segs$n.probes),
    end = rep(comb_res_segs$end.pos, times = comb_res_segs$n.probes),
    copynumbers = comb_res_bins[[s]],
    segments = rep(comb_res_segs[[s]], times = comb_res_segs$n.probes)
  )
  
  # Run squaremodel
  result_sqm_joint <- squaremodel(template_joint, cellularities = seq(95, 100, by = 0.1), ptop = 3, pbottom = 1.5, prows = 20)
  
  # Include number of reads
  used_reads <- number_of_reads$UsedReads[number_of_reads$Sample == s]
  
  # Get adjusted segments
  segmentdf_joint = getadjustedsegments(template_joint, cellularity = result_sqm_joint$minimadf[sol_number,"cellularity"], 
                                        ploidy = result_sqm_joint$minimadf[sol_number,"ploidy"])
  psit  = result_sqm_joint$minimadf[sol_number,"ploidy"]
  rho   = result_sqm_joint$minimadf[sol_number,"cellularity"]
  segmentdf_joint$Chromosome <- as.numeric(segmentdf_joint$Chromosome)
  segmentdf_joint <- segmentdf_joint %>% arrange(Chromosome, Start)
  
  # Plot
  # Prepare df for plotting
  plt.df_joint = data.frame(genome.bin = 1:nrow(na.omit(template_joint)), 
                            chr = factor(paste0("chr",rep(segmentdf_joint$Chromosome,
                                                          times = segmentdf_joint$Num_Bins)),
                                         levels = paste0("chr",c(1:22))), 
                            start = na.omit(template_joint)$start,
                            end = na.omit(template_joint)$end,
                            Log2ratio = log(na.omit(template_joint)$copynumbers, base = 2),
                            Call = as.character(rep(segmentdf_joint$Copies, times = segmentdf_joint$Num_Bins)),
                            Segment = log(na.omit(template_joint)$segments, base = 2))
  
  # Make the plot
  plot_joint = ggplot(plt.df_joint, aes(x = genome.bin, y = Log2ratio, col = Call)) +
    geom_hline(yintercept = c(-3,-2,-1,0,1,2,3), lty = c("solid"), lwd = 0.2) +
    geom_point() +
    ylab("log2ratio") +
    scale_colour_manual(values = cols[sort(unique(plt.df_joint$Call))]) +
    scale_x_continuous(name = "Chromosomes", labels = c(1:22), 
                       breaks = as.vector(c(1, cumsum(table(plt.df_joint$chr))[-22]) + 
                                            (table(plt.df_joint$chr) / 2))) + 
    geom_vline(xintercept = c(1, cumsum(table(plt.df_joint$chr))), lty = "dotted") +
    ggtitle(paste0(s, ", cellularity=",
                   result_sqm_joint$minimadf[sol_number,"cellularity"],", ploidy=",
                   result_sqm_joint$minimadf[sol_number,"ploidy"],", used reads=",used_reads)) + 
    scale_y_continuous(limits=c(-3,3), oob=scales::squish) + 
    geom_point(aes(y = Segment), color="#000000", size = 0.1) +
    #theme_cowplot() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14))
  
  # print(plot_joint)
  
  print(length(plt.df_joint$Call))
  
  #dir.create(paste0(out_dir, "/plots_joint"))
  ggsave(paste0(img_path,'/', s,"_", bin_size,"kbp_ACE_cn_profile_jointseg.png"),plot = plot_joint, height = 5, width = 15)
  
  plt.df_joint = plt.df_joint %>% 
    mutate(sample = s)
  
  return(list(plot = plot_joint, cna_calls = plt.df_joint))
  
})
names(joint_seg_res) = sample_names

joint_seg_plt = grid.arrange(joint_seg_res$A673_EFEndo_Rabbit_IgG_R1$plot, joint_seg_res$A673_EFEndo_Rabbit_IgG_R2$plot)
ggsave(paste0(img_path, '/EFEndo_IgG_joint_seg.png'), plot = joint_seg_plt, height = 8, width = 15, dpi = 300, units = 'in')

joint_seg_plt = grid.arrange(joint_seg_res$A673_EFKD_Rabbit_IgG_R1$plot, joint_seg_res$A673_EFKD_Rabbit_IgG_R2$plot)
ggsave(paste0(img_path, '/EFKD_IgG_joint_seg.png'), plot = joint_seg_plt, height = 8, width = 15, dpi = 300, units = 'in')

joint_seg_plt = grid.arrange(joint_seg_res$A673_EFRescue_Rabbit_IgG_R1$plot, joint_seg_res$A673_EFRescue_Rabbit_IgG_R2$plot)
ggsave(paste0(img_path, '/EFRescue_IgG_joint_seg.png'), plot = joint_seg_plt, height = 8, width = 15, dpi = 300, units = 'in')

joint_seg_plt = grid.arrange(joint_seg_res$A673_Rb_IgG_R1$plot, joint_seg_res$A673_Rb_IgG_R2$plot)
ggsave(paste0(img_path, '/EF_IgG_joint_seg.png'), plot = joint_seg_plt, height = 8, width = 15, dpi = 300, units = 'in')

