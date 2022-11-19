#' @export
#' @import methods
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
setMethod("seqinfo", signature("GPosExperiment"), 
          function(x) x@rowRanges@seqinfo)

#' Get  GPosExperiment scores vector
#' 
#' TODO Long Title
#' 
#' TODO Narraitve
#' 
#' @param x
#' @param apply_mask
#' @param zero_fill
#' 
#' @export
#' @import methods
setMethod("scores", 
  signature(x = "GPosExperiment"), 
  function(x, apply_mask = TRUE, zero_fill = TRUE) {
    stopifnot(isa(apply_mask, "logical"), 
              isa(zero_fill, "logical"))
    rows <- rowRanges(x)
    result <- lapply(colData(x)$NETseqData, function(u) {
      
      cols <- .process_colData_scores(u, 
                                      zero_fill = zero_fill, 
                                      apply_mask = apply_mask, 
                                      range_scope = rows)
  
      v <- vector("list", length(rows))
      ov <- findOverlaps(rows, cols)
      y <- split(cols[subjectHits(ov)], queryHits(ov))
      v[as.integer(names(y))] <- as.list(y)
    })
    matrix(unlist(result, recursive = FALSE), nrow = nrow(x), ncol = ncol(x)
           , dimnames = dimnames(x)
           )
  }
)

#' 
#' @export
#' @import methods
#' @importMethodsFrom GenomicRanges findOverlaps
setMethod("mask", 
          signature("GPosExperiment"), 
  function(x) {
    z <- lapply(x@colData$NETseqData, function(u) {
      v <- rep(list(GRanges()), length = nrow(x))
      ov <- findOverlaps(rowRanges(x), u@mask)
      ov <- lapply(split(subjectHits(ov), queryHits(ov)), 
                   function(w) list(u@mask[w])[[1]])
      v[as.integer(names(ov))] <-  ov
      v
    })
    matrix(unlist(z, recursive = FALSE), nrow = nrow(x), ncol = ncol(x)
           , dimnames = dimnames(x))
  })

#'
#' @export
#' @import methods
setMethod("plot", signature("GPosExperiment"),
          function(x) {
            .plot(x)
          })
