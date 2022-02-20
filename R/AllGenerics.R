#' @export
setGeneric("scores", function(x) standardGeneric("scores"))

#' @export
setGeneric("scores<-", function(x, value) standardGeneric("scores<-"))

#' @export
setGeneric("subranges", function(x) standardGeneric("subranges"))

#' @export
setGeneric("subranges<-", function(x, value) standardGeneric("subranges<-"))

#' @export
setGeneric("NETseqDataFromBedgraph",
       function(sampleId, ...)
       standardGeneric("NETseqDataFromBedgraph"))

#' @export
setGeneric("NETseqDataFromBAM",
           function(sampleId, ...)
             standardGeneric("NETseqDataFromBAM"))
