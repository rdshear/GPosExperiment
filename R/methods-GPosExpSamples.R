### =========================================================================
### GPosExpSamples objects
### -------------------------------------------------------------------------
###


# TODO: Create accessor functions for these slots (maybe indexed in parent class?)

setValidity("GPosExpSamples", function(object)
{
  msg <- NULL

  if (length(object@sampleNames) != length(object@scoreFilesPos) |
      length(object@sampleNames) != length(object@scoreFilesNeg)) {
    msg <- c(msg, gettextf("There must be the same number of score file names as there are sample names"))
  } else if (length(object@sampleNames) != length(object@segmentFiles)) {
      msg <- c(msg, gettextf("There must be the same number of score file names as there are sample names"))
    } else {
    for (f in object@scoreFilesPos) {
      x <- file.path(object@scoreFileDirectory, f)
      if (!file.exists(x)) {
        msg <- c(msg,
                 gettextf("(+) strand score file not found '%s'", x))
      }
    }

      for (f in object@scoreFilesNeg) {
        x <- file.path(object@scoreFileDirectory, f)
        if (!file.exists(x)) {
          msg <- c(msg,
                   gettextf("(-) strand score file not found '%s'", x))
        }
      }

      for (f in object@segmentFiles) {
        x <- file.path(object@segmentFileDirectory, f)
        if (!file.exists(x)) {
          msg <- c(msg,
                   gettextf("segment file not found '%s'", x))
        }
      }
  }

  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})

# TODO: sampleNames is a bioBase generic function...pick another name
setMethod("initialize",
          signature(.Object = "GPosExpSamples"),
          function(.Object,
                   scoreFileDirectory = "",
                   segmentFileDirectory = "",
                   sampleNames = character(),
                   scoreFilesPos = character(),
                   scoreFilesNeg = character(),
                   segmentFiles = character(),
                    ...)
          {
            .Object@scoreFileDirectory <- scoreFileDirectory
            .Object@segmentFileDirectory <- segmentFileDirectory
            .Object@sampleNames <- sampleNames
            .Object@scoreFilesPos <- scoreFilesPos
            .Object@scoreFilesNeg <- scoreFilesNeg
            .Object@segmentFiles <- segmentFiles
            validObject(.Object)
            .Object
          }
)

#' @exportClass GPosExpSamples
#' @export
GPosExpSamples <- function(scoreFileDirectory = "",
                       segmentFileDirectory = "",
                       sampleNames = character(),
                       scoreFilesPos = character(),
                       scoreFilesNeg = character(),
                       segmentFiles = character(),
                        ...) {
  new("GPosExpSamples",
      scoreFileDirectory = scoreFileDirectory,
      segmentFileDirectory = segmentFileDirectory,
      sampleNames = sampleNames,
      scoreFilesPos = scoreFilesPos,
      scoreFilesNeg = scoreFilesNeg,
      segmentFiles = segmentFiles,
      ...)
}

