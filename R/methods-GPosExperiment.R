#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
GPosExperiment <- function(sample = GPosExpSamples(),
                   seqinfo =  GenomeInfoDb::Seqinfo(),
                   rowRanges = GRanges(),
                   occupancyAssayName = "occupancy",
                   segmentsAssayName = "segments",
                   ...)
{

  # TODO: Move to initialize method
  # TODO: build the colData with a list, one element per row
  columns <- DataFrame(sample = sample@sampleNames,
                       score_file_pos =
                         file.path(sample@scoreFileDirectory, sample@scoreFilesPos),
                       score_file_neg =
                         file.path(sample@scoreFileDirectory, sample@scoreFilesNeg),
                       segment_file =
                         file.path(sample@segmentFileDirectory, sample@segmentFiles),
                       row.names = sample@sampleNames)

  object <- .GPosExperiment(SummarizedExperiment(rowRanges = rowRanges, colData = columns, ...),
                    occupancyAssayName = occupancyAssayName,
                    segmentsAssayName = segmentsAssayName)

  validObject(object)
  object
}

S4Vectors::setValidity2("GPosExperiment", function(object) {
  msg <- NULL

  # TODO: add validity
  # if (assayNames(object)[1] != "counts") {
  #   msg <- c(msg, "'counts' must be first assay")
  # }
  #
  # if (min(assay(object)) < 0) {
  #   msg <- c(msg, "negative values in 'counts'")
  # }
  #
  if (is.null(msg)) {
    TRUE
  } else msg
})


#' @import GenomicRanges
#' @import BiocIO
#' @import SummarizedExperiment
#' @import IRanges
#' @import S4Vectors
#' @import Rsamtools
#' @export
setMethod("readOccupancy", signature("GPosExperiment"), function(x, varname = "score") {
# TODO: currently returns an GPosExperiment argument
# Fix it so that behaves like CAGEr::getCTSS

  if (!x@occupancyAssayName %in% names(x@assays)) {
    g <- rowRanges(x)
    z <- mapply(function(pos_file, neg_file) {
      v <- mapply(function(f, s, query) {
              w <- import(f, which = query, seqinfo = seqinfo(x))
              w <- trim(w)
              w <- unlist(mcolAsRleList(w, varname = varname))
              gp <- GPos(seqinfo(x))
              strand(gp) <- s
              gp$occupancy <- Rle((function(z) ifelse(is.na(z), 0L, as.integer(z)))(runValue(w)), runLength(w))
              gp
            },
            c(pos_file, neg_file),
            c('+', '-'),
            split(g, strand(g))[1L:2L],
            USE.NAMES = FALSE, SIMPLIFY = FALSE)

      v <- c(v[[1L]], v[[2L]])
      w <- findOverlaps(g, v)
      list(lapply(split(subjectHits(w), queryHits(w)), function(u) v[u]))
    },
    colData(x)$score_file_pos,
    colData(x)$score_file_neg,
    USE.NAMES = FALSE)
    assay(x, x@occupancyAssayName) <- matrix(do.call(c, z), nrow = nrow(x))

    # TODO: Be sure to document this as not thread safe or add locks
    # overwrite the object in the parent environment.
    # TODO: Is this an acceptable practice?
    # The CAGEr method assign(objectName, x, parent.env()) doesn't seem to work

  }

  x
})


#' @export
setMethod("occupancy", signature("GPosExperiment"), function(x, ...) assay(x, x@occupancyAssayName))


#' @export
setMethod("occupancyRle", signature("GPosExperiment"), function(x, reverse.negative = FALSE, ...) {
  m <- occupancy(x, ...)
  # TODO: typechecks
  if (reverse.negative) {
    f <- sapply(strand(rowRanges(x)) == "-", function(w) ifelse(w, function(u) rev(u), function(u) u))
  }  else
  {
    f <- list(function(u) u)
  }

  m1 <- matrix(mapply(function(a, b) a(b$occupancy), f, m), nrow = nrow(x), ncol = ncol(x))
  dimnames(m1) <- dimnames(m)
  m1
})

#' seqinfo
#'
#' TODO: describe
#'
#' @export
#'
setMethod("seqinfo", signature("GPosExperiment"), function(x) x@rowRanges@seqinfo)


#' segments
#'
#' TODO: describe
#'

#' @export
setMethod("segments", signature("GPosExperiment"), function(x, reverse.negative = FALSE, ...) {
# TODO: merge with occupancyRle code in private function
  m <- assay(x, "segments")
  if (reverse.negative) {
    f <- sapply(strand(rowRanges(x)) == "-", function(w) ifelse(w, function(u) rev(u), function(u) u))
  }  else
  {
    f <- list(function(u) u)
  }

  m1 <- matrix(mapply(function(a, b) a(b), f, m), nrow = nrow(x), ncol = ncol(x))
  dimnames(m1) <- dimnames(m)
  m1
})
#' @import SummarizedExperiment
#' @export
setMethod("readSegments", signature("GPosExperiment"), function(x, ...) {
  # TODO: See comment at occupancy
  if (!x@segmentsAssayName %in% names(x@assays)) {
    g <- rowRanges(x)
    z <- lapply(colData(x)$segment_file, function(f) {
      v <- import(f, which = g)
      # TODO: Edit for needed mcol's
      mcols(v) <- NULL
      w <- findOverlaps(g, v)
      lapply(split(subjectHits(w), queryHits(w)), function(u) v[u])
    })

    assay(x, x@segmentsAssayName) <- matrix(do.call(c, z), nrow = nrow(x), ncol = ncol(x))

    # TODO: Be sure to document this as not thread safe or add locks
    # overwrite the object in the parent environment.
    # TODO: Is this an acceptable practice?
    # The CAGEr method assign(objectName, x, parent.env()) doesn't seem to work
    }
  x
})



# setMethod("calccp", "GPosExperiment", function(.Object, algorithm, f, parameters, f2) {
#   # TODO: Test types
#   fp <- formals(f)
#   for (w in names(parameters)) fp[[w]] <- parameters[[w]]
#   formals(f) <- fp
#   v <- genescores(.Object)
#   startTime <- Sys.time()
#   results <- data.frame(t(sapply(v, function(u) {
#     changeponts <- auxdata <- errors <- NULL
#     b <- try({
#       seg <- f(u)
#       c(1, f2(seg), length(u) + 1)
#     }, silent = TRUE)
#     if (inherits(b, c("try-error", "character"))) {
#       errors <- b
#     } else {
#       # TODO Declare error if > kMax returned
#       if (length(b) > .Object@kMax) {
#         b <- b[1:.Object@kMax]
#       }
#     changeponts <- b
#     }
#     # TODO: auxdata
#     c(changeponts = IntegerList(changeponts), auxdata = list(auxdata), errors = list(errors))
#   })))
#
#
#   endTime <- Sys.time()
#   runTime <- difftime(endTime, startTime, units = "secs")
#   logs <- c(runTime = as.character(runTime),
#             startTime = as.character(startTime),
#             endTime = as.character(endTime))
#   y <- GPosExpSamples(algorithm = algorithm,
#                      parameters = parameters, logs = logs, results = results)
#   y
# })

#TODO: Select which result(s) to plot
#       "Strips" for multiple results
#       Rationalize the alternate parameters
#       save / restore parameters
# plot.GPosExperiment <- function(object,
#                        filter = names(features(object)),
#                        tracks = seq_along(object@result.list),
#                        alpha = 0.6, color = "red", label = "", xaxis = TRUE, yaxis = TRUE, ...) {
#   oldpar <- par(no.readonly = TRUE)
#   if (length(filter) == 0) {
#     filter <- names(features(object))
#   }
#
#   gscore <- genescores(object, filter)
#   for (gene in filter) {
#     y <- gscore[[gene]]
#     ntracks <- length(tracks)
#     par(cex = 0.6, tcl = -0.25, mgp = c(2,0.6, 0))
#     par(mar = c(0.2, 0, 0, 0), oma = c(3, 3, 0.5, 0.5))
#     par(mfrow = c(ntracks, 1))
#
#     for (track in tracks) {
#       # TODO: create and use accessor function
#       r <- object@result.list[[track]]
#       cps <- r@results$changeponts
#       cp <- cps[[gene]][[1]]
#       # TODO: if error, display error instead of supressing everything
#       if (length(cp) > 0) {
#         begin <- cp[1:(length(cp) - 1)]
#         fini <- cp[2:length(cp)] - 1
#         seg.mean <- sapply(seq_along(begin), function(i) mean(y[begin[i]:fini[i]]))
#
#         plot(NULL, xlim = c(0,length(y)), ylim = range(log(y + 1)),
#              xaxt = ifelse(xaxis & track == tracks[ntracks], "s", "n"),
#              yaxt = ifelse(yaxis, "s", "n"), ...)
#         # TODO: figure out where to render 'label'
#         mtext(paste0(gene, ":", r@algorithm), side = 3, line = -1, adj = 0.015, cex = 0.7, font = 2)
#         polygon(cbind(c(rep(0:length(y), each = 2), 1), log(1 + c(0, rep(y, each = 2), 0, 0))),
#                 density = -1, col = scales::alpha(colour = color, alpha = alpha), border = NA )
#         abline(v = cp, col = "blue")
#         segments(x0 = begin, y0 = log(seg.mean + 1), x1 = fini, col = "black", lwd = 1)
#
#
#         mtext("base pairs downstream from TSS", side = 1, outer = TRUE, cex = 0.7, line = 2)
#         mtext(expression("log2(occupancy + 1)"), side = 2, outer = TRUE, cex = 0.7, line = 2)
#       }
#     }
#   }
#   par(oldpar)
# }

