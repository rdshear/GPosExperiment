#' Get test data filenames
#'
TestDataFilenames <- function() {
  fn <- Vectorize(function(x, ext) {
    system.file("extdata",
                paste0(x, ext, collapse = ""),
                package = "GPosExperiment")
  })
  
  sampleList <- c("SRR12840066", "SRR12840067")
  lapply(sampleList, function(u) {
    v <- fn(u, c(".pos.bedgraph.bgz",
                 ".neg.bedgraph.bgz",
                 ".bam"))
    data.frame(
      sampleId = u,
      bedgraph_pos = v[1],
      bedgraph_neg = v[2],
      bam = v[3]
    )
  }) %>%
    do.call(rbind, .) -> samples
  list(genes = fn("genelist", ".gff3.bgz"), seqinfo = Seqinfo(genome = "sacCer3") , samples = samples)
}
