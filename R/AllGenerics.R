#' @export
setGeneric("scores", function(x, ...) standardGeneric("scores"))

#' @export
setGeneric("scores<-", function(x, value) standardGeneric("scores<-"))

#' @export
setGeneric("mask", function(x) standardGeneric("mask"))

#' @export
setGeneric("mask<-", function(x, value) standardGeneric("mask<-"))

#' @export
setGeneric("HRseqDataFromBedgraph",
       function(sampleId, ...)
       standardGeneric("HRseqDataFromBedgraph"))

#' @export
setGeneric("HRseqDataFromBAM",
           function(sampleId, ...)
             standardGeneric("HRseqDataFromBAM"))
