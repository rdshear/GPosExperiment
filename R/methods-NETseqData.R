### =========================================================================
### NETseqData objects
### -------------------------------------------------------------------------
###


#' @import methods
setValidity("NETseqData", function(object)
{
  msg <- NULL
  # TODO Fill it in
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

# TODO see import.bedGraph genome parameter for name or object
##############
# function(con, format, text, trackLine = TRUE,
#          genome = NA, colnames = NULL,
#          which = NULL, seqinfo = NULL, extraCols = character(),
#          sep = c("\t", ""), na.strings=character(0L))
# {
#   if (!missing(format))
#     checkArgFormat(con, format)
#   sep <- match.arg(sep)
#   stopifnot(is.character(na.strings), !anyNA(na.strings))
#   file <- con
#   m <- manager()
#   con <- queryForConnection(m, con, which)
#   on.exit(release(m, con))
#   if (attr(con, "usedWhich"))
#     which <- NULL
#   if (is(genome, "Seqinfo")) {
#     seqinfo <- genome
#     genome <- NA_character_
#   }
#   if (is.null(seqinfo))
#     seqinfo <- attr(con, "seqinfo")
############
#' @import methods
#' @importClassesFrom GenomicRanges GRanges GPos
setMethod("initialize",
          signature(.Object = "NETseqData"),
          function(.Object,
                   scores = GPos(stitch = FALSE),
                  segments = GRanges(),
                  sampleId = character(),
                  seqinfo = NULL,
                 ...) 
  {
            .Object@scores <- scores
            .Object@segments <- segments
            .Object@sampleId <- sampleId
            .Object@seqinfo <- seqinfo
            validObject(.Object)
            .Object
          })

#' NETseqData Constructor
#' @param scores 
#' @param segments 
#' @param sampleId 
#' @param seqinfo
#' @param ... 
#'
#' @exportClass NETseqData
#' @export
#' @importClassesFrom GenomicRanges GRanges GPos
NETseqData <- function(scores = GPos(stitch = FALSE),
                          segments = GRanges(),
                          sampleId = character(),
                          seqinfo = NULL,
                      ...) {
  if (is.null(seqinfo)) {
    seqinfo <- seqinfo(segments)
  }
  new("NETseqData",
      scores = GPos(scores, stitch = FALSE),
      segments = segments,
      sampleId = sampleId,
      seqinfo = seqinfo,
      ...)
}

#' @export
setMethod("scores", signature(x = "NETseqData"), function(x) x@scores)

#' @export
setMethod("scores<-", signature(x = "NETseqData"), function(x, value) 
{
  x@scores <- value
  x
})

#' @export
setMethod("segments", signature(x = "NETseqData"), function(x) x@segments)

#' @export
setMethod("segments<-", signature(x = "NETseqData"), function(x, value) 
{
  x@segments <- value
  x
})

#' @export
setMethod("seqinfo", signature(x = "NETseqData"), function(x) x@seqinfo)

#' @export
setMethod("names", signature(x = "NETseqData"), function(x) x@sampleId)

#' @export
setMethod("names<-", signature(x = "NETseqData"), function(x, value) 
{
    x@sampleId <- value
    x
})

#' @export
#' 
#' @importMethodsFrom BiocIO import
#' @importClassesFrom GenomicRanges GRanges GPos
#' @importClassesFrom GenomeInfoDb Seqinfo
# TODO Add filter GRanges
setMethod("NETseqDataFromBedgraph", signature = c("character", "character"), 
  function(sampleId, filenames, seqinfo) {
    if (!isa(seqinfo, "Seqinfo")) {
      stop("seqinfo must be of type Seqinfo")
    }
    if (!is.character(filenames) || length(filenames) != 2 ||
        all(sort(names(filenames)) != c("-", "+"))) {
      stop("filenames must be character vector with elements named '+' and '-'")
    }
    x <- mapply(function(strand_sym, infilename) {
      x <- import(infilename, seqinfo = seqinfo)
      strand(x) <- strand_sym
      x
    }, list('+', '-'),  filenames, SIMPLIFY = FALSE)
    
    x <- GRangesToGPos(sort(c(x[[1]], x[[2]])))
    xs <- as.integer(x$score)
    if (all(x$score == xs)) {
      x$score <- xs
    }
    NETseqData(scores = x, sampleId = sampleId, seqinfo = seqinfo)
  })

### Helper FUnctions

#' GRangesToGPos
#'
#'
#' @param x 
#' @import BiocGenerics
GRangesToGPos <- function(x) {
  GPos(x, score = rep(x$score, width(x)))
}

# TODO show
# TODO plot
# TODO print