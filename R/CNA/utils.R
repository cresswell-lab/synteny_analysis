run_ace = FALSE
bin_size = 1000

# cell lines expected ploidies 

exp_psi = tibble(
  cell_line = c('KNS42', 'SF9427', "hNSC", "SF9402", "SF7761", "NHA", "DIPG_IV", "SF8628", "DIPG_XIII","DIPG007" ),
  ploidy = c(3, 2, 2, 2, 2, 2,2,2,2,2)
)

# Set parameters
cna_mat_joint = NULL
psit_mat = NULL
rho_mat = NULL
sol_number = 1

extract_results = function(x, 
                           genome = 'hg38', 
                           gamma = 40, 
                           ptop, 
                           pbottom, 
                           number_of_reads, 
                           cols, 
                           id) {
  
  # Calculate logged copynumbers and add position  
  x = x %>% 
    mutate(log2_copynumber = log(abs(copynumbers), base = 2)) %>% 
    mutate(pos = (end + start)/2) %>% 
    # Get rid of all the columns we do not use
    dplyr::select(chr, pos, Sample, log2_copynumber) %>%
    filter(chr != "Y")
  
  # Group to get the log2value of each sample at each position of each chromosome
  df <- x %>%
    group_by(chr, pos, Sample) %>%
    dplyr::slice(1) %>%  
    ungroup() %>% 
    mutate(chr = factor(chr, levels = c(seq(1:22), 'X'))) %>% 
    arrange(chr, pos)
  
  sample_id = df$Sample %>% unique

  # change format so that each sample is one column
  df_wide <- df %>%
    pivot_wider(
      names_from = Sample,
      values_from = log2_copynumber
    )
  
  # exclude NAs
  df_wide <- na.omit(df_wide) %>% 
    as.data.frame()
  
  res_segs = pcf(df_wide, assembly = genome, gamma = gamma)
  
  # Data preparation
  # res_segs contains the mean logged copynumbers per segment --> need the data not logged
  res_segs[,7] <- 2^(res_segs[,7])
  
  # single_purity_wide contains the log2 ratio for each sample per bin --> also need it not logged for the function
  res_bins <- df_wide
  res_bins[, c(3:ncol(res_bins))] <- 2^(df_wide[, c(3:ncol(res_bins))])
  
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
    copynumbers = res_bins[[sample_id]],
    segments = rep(res_segs[[sample_id]], times = res_segs$n.probes)
  )
  
  # Run squaremodel
  result_sqm <- squaremodel(template, cellularities = seq(95, 100, by = 0.1), ptop = ptop, pbottom = pbottom, prows = 20)
  #ploidy 1.5-2.5 for diploid, 2.5-3.5 for triploid and 3.5-4.5 for tetraploid
  
  # Include number of reads
  used_reads = number_of_reads %>% 
    filter(Sample == sample_id) %>% 
    pull(UsedReads)
  
  # Get adjusted segments
  segmentdf = getadjustedsegments(template, cellularity = result_sqm$minimadf[sol_number,"cellularity"], 
                                  ploidy = result_sqm$minimadf[sol_number,"ploidy"])
  psit  = result_sqm$minimadf[sol_number,"ploidy"]
  rho   = result_sqm$minimadf[sol_number,"cellularity"]
  segmentdf$Chromosome <- as.numeric(segmentdf$Chromosome)
  segmentdf <- segmentdf %>% arrange(Chromosome, Start)
  
  # Plot

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
    ggtitle(id) + 
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
  
  # dir.create(file.path(out_dir, "img"))
  # ggsave(file.path(out_dir, "img", paste0(sample, '_', bin_size, '_kbp_ACE_cn_profile.png')), 
  #        plot = plot, height = 4, width = 15, dpi = 300, units = 'in')
  #ggsave(paste0(out_dir,"/plots_thesis/",sample,"_",
  #bin_size,"kbp_ACE_cn_profile.png"), 
  #plot = plot,
  #height = 3.5, width = 15)
  
  # Collect calls
  plt.df$Call <- as.character(rep(segmentdf$Copies, times = segmentdf$Num_Bins)) 
  
  # cna_mat = cbind(cna_mat, plt.df$Call) 
  model_check = cbind(result_sqm)
  
  res = list(
    'squaremodel' = result_sqm,
    'psit_mat' = psit,
    'rho_mat' = rho, 
    'model_check' = model_check, 
    'df' = plt.df, 
    'plot' = plot
  )
  return(res)
  
  # saveRDS(file = paste0(out_dir, "/RDS_thesis/STA-NB-13/", sample, "_segments.rds"), psit_mat)
}

