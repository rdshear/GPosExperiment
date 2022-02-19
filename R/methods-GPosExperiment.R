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
setMethod("scores", signature(x = "GPosExperiment"), function(x) {
  r <- rowRanges(x)
  result <- lapply(colData(x)$NETseqData, function(u) {
    cols <- u@scores
    v <- vector("list", length(r))
    ov <- findOverlaps(rowRanges(x), cols)
    y <- lapply(split(ov, queryHits(ov)), function(u) 
      as.integer(cols[subjectHits(u)]$score))
    v[as.integer(names(y))] <- y
    v
  })
  matrix(unlist(result, recursive = FALSE), nrow = nrow(x), ncol = ncol(x))
})



