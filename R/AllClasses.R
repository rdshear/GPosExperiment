#' An S4 sublclass of RangedSummairzedExperiment
#' 
#' TBD
#' 
#' @slot seqinfo 
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment

.GPosExperiment <- setClass("GPosExperiment",
                    slots = representation(
                      seqinfo = "Seqinfo"
                    ),
                    
                    contains = "RangedSummarizedExperiment"
)

#' Contains sample-level assay date
#' 
#' TBD
#'
#' @slot sampleId 
#' @slot scores 
#' @slot subranges 
#' @slot seqinfo 
#'
#' @import methods
#' @importClassesFrom GenomicRanges GRanges GPos
#' @importClassesFrom GenomeInfoDb Seqinfo
setClass("NETseqData",
                        slots = c(
                          sampleId = "character",
                          scores = "GPos",
                          subranges = "GRanges",
                          seqinfo = "Seqinfo"
                        ))