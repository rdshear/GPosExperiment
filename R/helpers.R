
.complementStrand <- function(u) c(`+`="-", `-`="+")[as.character(strand(u))]

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
    GPos(., score = rep(.$V1, width(.)), seqinfo = seqinfo(u))
}

.OverlappedRanges <- function(q, s) {
  result <- s[subjectHits(findOverlaps(q, s))]
  seqinfo(result) <- seqinfo(s)
  result
}
# TODO Is this used
.GRangesToZeroFillGPos <- function(u) {
  # This function needs at least one element in u to function correctly
  # So if no elements in u, then create a single zero score so it can run
  if (length(u) == 0) {
    u <- append(u, 
                GRanges(seqnames = seqnames(seqinfo(u))[1], 
                  ranges = IRanges(start = 1), strand = "+", 
                  score = 0L))
  } 
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
  score_int <- as.integer(v$score)
  if (all(score_int == v$score))
  {
    v$score = score_int
  }
  GPos(v, score = rep(v$score, width(v)))
}

#' Get test data filenames
#'
.TestDataFilenames <- function() {
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
