#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom S4Vectors DataFrame
GPosExperiment <- function(sample = NETseqData(),
                   seqinfo =  GenomeInfoDb::Seqinfo(),
                   rowRanges = GRanges(),
                   ...)
{
  cols <- DataFrame(NETseqData = List(sample), row.names = names(sample))
  
  object <- .GPosExperiment(SummarizedExperiment(rowRanges = rowRanges, 
                       colData = cols, ...))

  validObject(object)
  object
}

#' @import methods
setValidity("GPosExperiment", function(object) {
  msg <- NULL

  # TODO: add validity
  if (is.null(msg)) {
    TRUE
  } else msg
})

#' seqinfo
#'
#' TODO: describe
#'
#' @export
#'
setMethod("seqinfo", signature("GPosExperiment"), 
          function(x) x@rowRanges@seqinfo)

#' Get  GPosExperiment scores vector
#' 
#' TODO Long Title
#' 
#' TODO Narraitve
#' 
#' @param x
#' @param apply_mask
#' @param zero_fill
#' 
#' @export
#' @import methods
setMethod("scores", 
  signature(x = "GPosExperiment"), 
  function(x, apply_mask = TRUE, zero_fill = TRUE) {
    stopifnot(isa(apply_mask, "logical"), 
              isa(zero_fill, "logical"))
    rows <- rowRanges(x)
    result <- lapply(colData(x)$NETseqData, function(u) {
      
      cols <- .process_colData_scores(u, 
                                      zero_fill = zero_fill, 
                                      apply_mask = apply_mask, 
                                      range_scope = rows)
  
      v <- vector("list", length(rows))
      ov <- findOverlaps(rows, cols)
      y <- split(cols[subjectHits(ov)], queryHits(ov))
      v[as.integer(names(y))] <- as.list(y)
    })
    matrix(unlist(result, recursive = FALSE), nrow = nrow(x), ncol = ncol(x)
           , dimnames = dimnames(x)
           )
  }
)

#' 
#' @export
#' @import methods
#' @importMethodsFrom GenomicRanges findOverlaps
setMethod("mask", 
          signature("GPosExperiment"), 
  function(x) {
    z <- lapply(x@colData$NETseqData, function(u) {
      v <- rep(list(GRanges()), length = nrow(x))
      ov <- findOverlaps(rowRanges(x), u@mask)
      ov <- lapply(split(subjectHits(ov), queryHits(ov)), 
                   function(w) list(u@mask[w])[[1]])
      v[as.integer(names(ov))] <-  ov
      v
    })
    matrix(unlist(z, recursive = FALSE), nrow = nrow(x), ncol = ncol(x)
           , dimnames = dimnames(x))
  })

#TODO: Select which result(s) to plot
#       "Strips" for multiple results
#       Rationalize the alternate parameters
#       save / restore parameters
# plot.HRESseq <- function(object,
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

