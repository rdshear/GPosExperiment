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
  # if (assayNames(object)[1] != "counts") {
  #   msg <- c(msg, "'counts' must be first assay")
  # }
  #
  # if (min(assay(object)) < 0) {
  #   msg <- c(msg, "negative values in 'counts'")
  # }
  #
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
  col_scores <- apply(x@colData, 1, function(u) (u$NETseqData)@scores)
  cols <- col_scores[[1]]
  ov <- findOverlaps(rowRanges(x), cols)
  # TODO THIS DOES NOt WORK if there are no overlaps for a specific cell
  # ...ransform to [[x,y]] addressing format (keep the single findOverlaps for performance)
  y <- lapply(split(ov, queryHits(ov)), function(u) 
    as.integer(cols[subjectHits(u)]$score))
  matrix(y, ncol = ncol(x), nrow = nrow(x))
})



