
#' @import methods
setValidity("NETseqData", function(object)
{
  msg <- NULL
  if (length(object@sampleId) != 1 || object@sampleId == "")
  {
    msg <- append(msg, "sampleId missing")
  }
  
  if (!isa(object@scores, "GPos")) {
    msg <- append(msg, "scores object must be a GPos object")
  }
  
  # TODO Reinstate?
  # if (length(seqinfo(object@scores)) < 1) {
  #   msg <- append(msg, "scores object must have a non-empty Seqinfo object")
  # }
  
  if (!isa(object@mask, "GRanges")) {
    msg <- append(msg, "mask object must be a GRanges object")
  }

  if (length(object@mask) > 0) {
    if (length(seqinfo(object@mask)) < 1) {
      msg <- append(msg, "mask object must have a non-empty Seqinfo object")
    }
  
    if (!identical(seqinfo(object@scores), seqinfo(object@mask))) {
      msg <- append(msg, "scores and mask must have the same seqinfo object")
    }
  } 
  
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})

#' Optionally Apply masks to a scores GPos and/or zero fill the gaps
#' 
.process_colData_scores <- function(u, apply_mask, zero_fill, range_scope = NULL) {
  y <- u@scores

  if (is.null(range_scope)) {
    q = unique(DataFrame(chr = as.character(seqnames(y)), s = strand(y)))
    range_scope <- unlist(GRangesList(lapply(split(q, q$s), function(v) {
      result <- GRanges((seqinfo(u)[v$chr]))
      strand(result) <- v$s
      result
    })))
  }
  y <- .OverlappedRanges(range_scope, y)
  
  if (zero_fill) {
    g <- GenomicRanges::intersect(range_scope, gaps(GRanges(y)))
    g <- g[as.character(strand(g)) != "*"]
    g$score <- 0L
    y <- c(GRangesToGPos(g), y)
  }

  if (apply_mask) {
    ov <- findOverlaps(u@mask, y)
    if (length(ov) > 0) {
      y[subjectHits(ov)]$score <- NA
    }
  }
  sort(y)
}

#' @import methods
#' @importClassesFrom GenomicRanges GRanges GPos
setMethod("initialize",
          signature(.Object = "NETseqData"),
          function(.Object,
                   scores = GPos(),
                  mask = GRanges(),
                  sampleId = character(),
                  seqinfo = NULL,
                
                 ...)
  {
            if (isa(scores, "GRanges") && !isa(scores, "GPos")) {
              scores <- GRangesToGPos(scores)
            }
            .Object@scores <- scores
            .Object@mask <- mask
            .Object@sampleId <- sampleId
            .Object@seqinfo <- seqinfo
            validObject(.Object)
            .Object
          })

#' NETseqData - Construct NETseqData
#' 
#' TBD
#' 
#' @param scores 
#' @param mask 
#' @param sampleId 
#' @param seqinfo
#' @param ... 
#'
#' @exportClass NETseqData
#' @export
#' @importClassesFrom GenomicRanges GRanges GPos
NETseqData <- function(scores = GRanges(),
                          mask = GRanges(),
                          sampleId = character(),
                          seqinfo = NULL,
                      ...) {
  if (!isa(scores, "GRanges") && !isa(scores, "GPos")) {
    stop("scores parameters must be GRanges or GPos")
  }  
  if (is.null(seqinfo)) {
    seqinfo <- seqinfo(scores)
  }
  if (is.null(seqinfo(scores))) {
    seqinfo(scores) <- seqinfo
  }
  new("NETseqData",
      scores = scores,
      mask = mask,
      sampleId = sampleId,
      seqinfo = seqinfo,
      ...)
}

#' scores -   NETseqData scores vector
#' 
#' TODO Long Title
#' 
#' TODO Narraitve
#' 
#' @param x
#' @param apply_mask
#' @param zero_fill
#' @param range_scope
#' 
#' @export
setMethod("scores", 
  signature(x = "NETseqData"), 
  function(x, apply_mask = FALSE, 
           zero_fill = FALSE, 
           range_scope = NULL) {
    stopifnot(isa(apply_mask, "logical"), 
              isa(zero_fill, "logical"),
              is.null(range_scope) || isa(range_scope, "GRanges"))
  
    .process_colData_scores(x, apply_mask, zero_fill, range_scope)
  }
)

#' @export
setMethod("scores<-", signature(x = "NETseqData"), function(x, value) 
{
  x@scores <- value
  x
})

#' @export
setMethod("mask", signature(x = "NETseqData"), function(x) x@mask)

#' @export
setMethod("mask<-", signature(x = "NETseqData"), function(x, value) 
{
  x@mask <- value
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
      NETseqData(scores = x, sampleId = s, seqinfo = seqinfo, mask = y)
    }, sampleId, filename_pos, filename_neg, filename_seg, SIMPLIFY = FALSE)
    names(result) <- sampleId
    result
  })

### Helper Functions

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