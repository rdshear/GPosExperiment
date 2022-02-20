#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment

.GPosExperiment <- setClass("GPosExperiment",
                    slots = representation(
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
                          scores = "GPos",
                          subranges = "GRanges",
                          seqinfo = "Seqinfo"
                        ))