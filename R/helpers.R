
complementStrand <- function(u) c(`+`="-", `-`="+")[as.character(strand(u))]

SamToScore <- function(u) {
  split(u, strand(u))[1:2] %>% as.list %>%
    map2(names(.), function(u, s) {
      coverage(u) %>% 
        bindAsGRanges %>%
        {strand(.) <- s; .}
    }) %>%
    GRangesList %>% 
    unlist %>%
    sort.GenomicRanges %>%
    GPos(., score = rep(.$V1, width(.)))
}

OverlappedRanges <- function(q, s) s[subjectHits(findOverlaps(q, s))]

GRangesToZeroFillGPos <- function(u) {
  w <- lapply(split(GRanges(u), strand(u)), mcolAsRleList, "score")
  v <- sapply(w, function(q) {
    y <- lapply(q, function(z) {
      runValue(z)[is.na(runValue(z))] <- 0
      z
    })
    y <- bindAsGRanges(score = as(y, "RleList"))
  })[c("+","-")]
  for (i in seq_along(v)) {
    strand(v[[i]]) <- names(v)[i]
  }
  v <- c(v[[1]], v[[2]])
  GPos(v, score = rep(v$score, width(v)))
}

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
