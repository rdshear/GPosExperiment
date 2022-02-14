#' @import methods
setClass("GPosExpSamples",
         slots = c(scoreFileDirectory = "character", # the file path for the samples
                   segmentFileDirectory = "character", # directory path for the changepoit files (gff3 format)
                   sampleNames = "character", # A vector of sample names
                   scoreFilesPos = "character", # Vector of filenames of + strand score files (bedGraph)
                   scoreFilesNeg = "character", # Vector of filenames of - strand score files (bedGraph)
                   segmentFiles = "character" #Vector of filenames for the segment files
         )
)


#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment

.GPosExperiment <- setClass("GPosExperiment",
                    slots = representation(
                      occupancyAssayName = "character",
                      segmentsAssayName = "character",
                      sampleData = "GPosExpSamples"
                    ),
                    
                    contains = "RangedSummarizedExperiment"
)
