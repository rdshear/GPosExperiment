#' @export
setGeneric("scores", function(x) standardGeneric("scores"))

#' @export
setGeneric("vscores", function(x, ...) standardGeneric("vscores"))

#' @export
setGeneric("scores<-", function(x, value) standardGeneric("scores<-"))

#' @export
setGeneric("mask", function(x) standardGeneric("mask"))

#' @export
setGeneric("mask<-", function(x, value) standardGeneric("mask<-"))

#' @export
setGeneric("NETseqDataFromBedgraph",
       function(sampleId, ...)
       standardGeneric("NETseqDataFromBedgraph"))

#' @export
setGeneric("NETseqDataFromBAM",
           function(sampleId, ...)
             standardGeneric("NETseqDataFromBAM"))
