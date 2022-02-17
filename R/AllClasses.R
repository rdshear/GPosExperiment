#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment

.GPosExperiment <- setClass("GPosExperiment",
                    slots = representation(
                      occupancyAssayName = "character",
                      segmentsAssayName = "character",
                      seqinfo = "Seqinfo"
                    ),
                    
                    contains = "RangedSummarizedExperiment"
)

#' @import methods
#' @importClassesFrom GenomicRanges GRanges GPos
#' @importClassesFrom GenomeInfoDb Seqinfo
setClass("NETseqData",
                        slots = c(
                          sampleId = "character",
                          scores = "UnstitchedGPos",
                          segments = "GRanges",
                          seqinfo = "Seqinfo"
                        ))