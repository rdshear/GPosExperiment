### =========================================================================
### NETseqData objects
### -------------------------------------------------------------------------
###


#' @import methods
setValidity("NETseqData", function(object)
{
  msg <- NULL
  # TODO Fill it in
  # todo verify that genelist is named

  
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
                   scores = GPos(stitch = FALSE),
                  subranges = GRanges(),
                  sampleId = character(),
                  seqinfo = NULL,
                
                 ...) 
  {
            .Object@scores <- scores
            .Object@subranges <- subranges
            .Object@sampleId <- sampleId
            .Object@seqinfo <- seqinfo
            validObject(.Object)
            .Object
          })

#' NETseqData Constructor
#' @param scores 
#' @param subranges 
#' @param sampleId 
#' @param seqinfo
#' @param ... 
#'
#' @exportClass NETseqData
#' @export
#' @importClassesFrom GenomicRanges GRanges GPos
NETseqData <- function(scores = GPos(stitch = FALSE),
                          subranges = GRanges(),
                          sampleId = character(),
                          seqinfo = NULL,
                      ...) {
  if (is.null(seqinfo)) {
    seqinfo <- seqinfo(subranges)
  }
  if (isa(scores, "GRanges")) {
    scores <- GRangesToGPos(scores)
  } else
  if (!isa(scores, "GPos")) {
    stop("scores parameters must be GRanges or GPos")
  }  
  new("NETseqData",
      scores = GPos(scores, stitch = FALSE),
      subranges = subranges,
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
setMethod("subranges", signature(x = "NETseqData"), function(x) x@subranges)

#' @export
setMethod("subranges<-", signature(x = "NETseqData"), function(x, value) 
{
  x@subranges <- value
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
#' @importMethodsFrom rtracklayer import
#' @importClassesFrom GenomicRanges GRanges GPos
#' @importClassesFrom GenomeInfoDb Seqinfo
# TODO Add filter GRanges
setMethod("NETseqDataFromBedgraph", signature = c("character"), 
  function(sampleId, filename_pos, filename_neg, filename_seg = NA, seqinfo) {
    if (length(sampleId) != length(filename_pos) ||
        length(sampleId) != length(filename_neg)) {
      stop("sampleId, filename_pos, filename_neg must all have same length")
    }
    if (!is.na(filename_seg) &&  length(sampleId) != length(filename_seg)) {
      stop("filename_seg must if same length of sampleId if present")
    }
    if (!isa(seqinfo, "Seqinfo")) {
      stop("seqinfo must be of type Seqinfo")
    }
    result <- mapply(function(s, fp, fn, seg) {
      x <- mapply(function(strand_sym, infilename) {
        x <- import(infilename, seqinfo = seqinfo)
        strand(x) <- strand_sym
        x
      }, list('+', '-'),  c(fp, fn), SIMPLIFY = FALSE)
      
      x <- GRangesToGPos(sort(c(x[[1]], x[[2]])))
      xs <- as.integer(x$score)
      if (all(x$score == xs)) {
        x$score <- xs
      }
      # TODO seqinfo -> genome names?
      if  (is.na(seg)) {
        y <-  GRanges()
      } else {
          y <- import(seg)
      }
      seqinfo(y) <- seqinfo
      NETseqData(scores = x, sampleId = s, seqinfo = seqinfo, subranges = y)
    }, sampleId, filename_pos, filename_neg, filename_seg, SIMPLIFY = FALSE)
    names(result) <- sampleId
    result
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