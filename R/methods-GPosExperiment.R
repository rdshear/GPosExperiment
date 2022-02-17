#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
GPosExperiment <- function(sample = NETseqData(),
                   seqinfo =  GenomeInfoDb::Seqinfo(),
                   rowRanges = GRanges(),
                   occupancyAssayName = "occupancy",
                   segmentsAssayName = "segments",
                   ...)
{
  cols <- DataFrame(NETseqData = List(sample), row.names = names(sample))
  
  object <- .GPosExperiment(SummarizedExperiment(rowRanges = rowRanges, colData = cols, ...),
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

