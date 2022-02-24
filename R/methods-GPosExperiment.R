#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
GPosExperiment <- function(sample = NETseqData(),
                   seqinfo =  GenomeInfoDb::Seqinfo(),
                   rowRanges = GRanges(),
                   ...)
{
  cols <- DataFrame(NETseqData = List(sample), row.names = names(sample))
  
  object <- .GPosExperiment(SummarizedExperiment(rowRanges = rowRanges, 
                       colData = cols, ...))

  validObject(object)
  object
}

#' @import methods
setValidity("GPosExperiment", function(object) {
  msg <- NULL

  # TODO: add validity
  if (is.null(msg)) {
    TRUE
  } else msg
})

#' seqinfo
#'
#' TODO: describe
#'
#' @export
#'
setMethod("seqinfo", signature("GPosExperiment"), function(x) x@rowRanges@seqinfo)

#' @export
setMethod("vscores", 
          signature(x = "GPosExperiment"), 
          function(x, apply.mask = TRUE) {
  r <- rowRanges(x)
  # TODO Performance Problem...refactor
  result <- lapply(colData(x)$NETseqData, function(u) {
    cols <- u@scores
    mask_ranges <- u@mask
    rows <- rowRanges(x)
    v <- vector("list", length(r))
    ov <- findOverlaps(rows, cols)
    y <- lapply(split(ov, queryHits(ov)), function(v) {
      g <- rows[queryHits(v)[1]]
      s <- cols[subjectHits(v)]
      dj <- disjoin(c(g,s), with.revmap = TRUE)
      dj <- dj[sapply(dj$revmap, length) == 1]
      sz <- GPos(dj, score = rep(0L, sum(width(dj))))
      s <- sort(c(sz, s))
      if (apply.mask) {
        ov <- findOverlaps(u@mask, s)
        if (length(ov) > 0) {
          s[subjectHits(ov)]$score <- NA
        }      
      }
      s$score
    })
    v[as.integer(names(y))] <- y
    v
  })
  matrix(unlist(result, recursive = FALSE), nrow = nrow(x), ncol = ncol(x))
})



