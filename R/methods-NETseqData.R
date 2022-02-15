### =========================================================================
### NETseqData objects
### -------------------------------------------------------------------------
###


# TODO: Create accessor functions for these slots (maybe indexed in parent class?)
#' @import methods
setValidity("NETseqData", function(object)
{
  msg <- NULL
  
  # if (length(object@sampleNames) != length(object@scoreFilesPos) |
  #     length(object@sampleNames) != length(object@scoreFilesNeg)) {
  #   msg <- c(msg, gettextf("There must be the same number of score file names as there are sample names"))
  # } else if (length(object@sampleNames) != length(object@segmentFiles)) {
  #   msg <- c(msg, gettextf("There must be the same number of score file names as there are sample names"))
  # } else {
  #   for (f in object@scoreFilesPos) {
  #     x <- file.path(object@scoreFileDirectory, f)
  #     if (!file.exists(x)) {
  #       msg <- c(msg,
  #                gettextf("(+) strand score file not found '%s'", x))
  #     }
  #   }
  #   
  #   for (f in object@scoreFilesNeg) {
  #     x <- file.path(object@scoreFileDirectory, f)
  #     if (!file.exists(x)) {
  #       msg <- c(msg,
  #                gettextf("(-) strand score file not found '%s'", x))
  #     }
  #   }
  #   
  #   for (f in object@segmentFiles) {
  #     x <- file.path(object@segmentFileDirectory, f)
  #     if (!file.exists(x)) {
  #       msg <- c(msg,
  #                gettextf("segment file not found '%s'", x))
  #     }
  #   }
  # }
  
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})

#' @import methods
#' @importClassesFrom GenomicRanges GRanges GPos
setMethod("initialize",
          signature(.Object = "NETseqData"),
          function(.Object,
                   scores = GRanges(),
                  segments = GRanges(),
                  sampleId = character(),
                 ...) 
  {
            .Object@scores <- scores
            .Object@segments <- segments
            .Object@sampleId <- sampleId
            validObject(.Object)
            .Object
          })

#' @param scores 
#'
#' @param segments 
#' @param sampleId 
#' @param ... 
#'
#' @exportClass NETseqData
#' @export
#' @importClassesFrom GenomicRanges GRanges GPos
NETseqData <- function(scores = GPos(stitch = FALSE),
                        segments = GRanges(),
                        sampleId = character(),
                      ...) {
  new("NETseqData",
      scores = GPos(scores, stitch = FALSE),
      segments = segments,
      sampleId = sampleId,
      ...)
}
