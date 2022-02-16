#' @export
setGeneric("readOccupancy", function(x, ...) standardGeneric("readOccupancy"))

#' @export
setGeneric("occupancy", function(x, ...) standardGeneric("occupancy"))

#' @export
setGeneric("occupancyRle", function(x, ...) standardGeneric("occupancyRle"))

#' @export
setGeneric("segments", function(x, ...) standardGeneric("segments"))

#' @export
setGeneric("segments<-", function(x, value) standardGeneric("segments<-"))

#' @export
setGeneric("scores", function(x) standardGeneric("scores"))

#' @export
setGeneric("scores<-", function(x, value) standardGeneric("scores<-"))

#' @export
setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

#' @export
setGeneric("readSegments", function(x, ...) standardGeneric("readSegments"))

#' @export
setGeneric("NETseqDataFromBedgraph",
       function(sampleId, filenames, seqinfo)
       standardGeneric("NETseqDataFromBedgraph"))


#TODO Add value
# #TODO: parameter for isolation mode + parameter for parallel operation
# setGeneric("calccp", function(.Object, ...) 0)