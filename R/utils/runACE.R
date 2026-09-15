runACE = function (inputdir = "./", outputdir, filetype = "rds", genome = "hg19", 
                   binsizes, ploidies = 2, imagetype = "pdf", method = "RMSE", 
                   penalty = 0, cap = 12, bottom = 0, trncname = FALSE, printsummaries = TRUE, 
                   savereadcounts = FALSE, autopick = FALSE, isDup = FALSE) 
{
  imagefunction <- get(imagetype)
  if (substr(inputdir, nchar(inputdir), nchar(inputdir)) != 
      "/") {
    inputdir <- paste0(inputdir, "/")
  }
  if (missing(outputdir)) {
    outputdir <- substr(inputdir, 0, nchar(inputdir) - 1)
  }
  if (!dir.exists(outputdir)) {
    dir.create(outputdir)
  }
  if (filetype == "bam") {
    if (missing(binsizes)) {
      binsizes <- c(100, 500, 1000)
    }
    parameters <- data.frame(options = c("inputdir", "outputdir", 
                                         "filetype", "binsizes", "ploidies", "imagetype", 
                                         "method", "penalty", "cap", "bottom", "trncname", 
                                         "printsummaries", "autopick"), values = c(inputdir, 
                                                                                   outputdir, filetype, paste0(binsizes, collapse = ", "), 
                                                                                   paste0(ploidies, collapse = ", "), imagetype, method, 
                                                                                   penalty, cap, bottom, trncname, printsummaries, 
                                                                                   autopick))
    for (b in binsizes) {
      currentdir <- file.path(outputdir, paste0(b, "kbp"))
      dir.create(currentdir)
      bins <- QDNAseq::getBinAnnotations(binSize = b, 
                                         genome = genome)
      readCounts <- QDNAseq::binReadCounts(bins, path = inputdir, isDuplicate = isDup)
      if (savereadcounts == TRUE) {
        saveRDS(readCounts, file = file.path(outputdir, 
                                             paste0(b, "kbp-raw.rds")))
      }
      readCountsFiltered <- QDNAseq::applyFilters(readCounts, 
                                                  residual = TRUE, blacklist = TRUE)
      readCountsFiltered <- QDNAseq::estimateCorrection(readCountsFiltered)
      copyNumbers <- QDNAseq::correctBins(readCountsFiltered)
      copyNumbers <- QDNAseq::normalizeBins(copyNumbers)
      copyNumbers <- QDNAseq::smoothOutlierBins(copyNumbers)
      copyNumbersSegmented <- QDNAseq::segmentBins(copyNumbers, 
                                                   transformFun = "sqrt")
      copyNumbersSegmented <- QDNAseq::normalizeSegmentedBins(copyNumbersSegmented)
      saveRDS(copyNumbersSegmented, file = file.path(outputdir, 
                                                     paste0(b, "kbp.rds")))
      ploidyplotloop(copyNumbersSegmented, currentdir, 
                     ploidies, imagetype, method, penalty, cap, bottom, 
                     trncname, printsummaries, autopick)
    }
  }
  else if (filetype == "rds") {
    parameters <- data.frame(options = c("inputdir", "outputdir", 
                                         "filetype", "ploidies", "imagetype", "method", "penalty", 
                                         "cap", "bottom", "trncname", "printsummaries", "autopick"), 
                             values = c(inputdir, outputdir, filetype, paste0(ploidies, 
                                                                              collapse = ", "), imagetype, method, penalty, 
                                        cap, bottom, trncname, printsummaries, autopick))
    files <- list.files(inputdir, pattern = "\\.rds$")
    for (f in seq_along(files)) {
      currentdir <- file.path(outputdir, paste0(substr(files[f], 
                                                       0, nchar(files[f]) - 4)))
      dir.create(currentdir)
      copyNumbersSegmented <- readRDS(file.path(inputdir, 
                                                files[f]))
      ploidyplotloop(copyNumbersSegmented, currentdir, 
                     ploidies, imagetype, method, penalty, cap, bottom, 
                     trncname, printsummaries, autopick)
    }
  }
  else {
    print("not a valid filetype")
  }
  write.table(parameters, file = file.path(outputdir, "parameters.tsv"), 
              quote = FALSE, sep = "\t", na = "", row.names = FALSE)
}
