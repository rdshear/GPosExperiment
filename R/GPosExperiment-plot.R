#TODO: Select which result(s) to plot
#       "Strips" for multiple results
#       Rationalize the alternate parameters
#       save / restore parameters

.plot <- function(object)
                       # filter = names(features(object)),
                       # tracks = seq_along(object@result.list),
                       # alpha = 0.6, color = "red", label = "", xaxis = TRUE, yaxis = TRUE, ...) 
{
  oldpar <- par(no.readonly = TRUE)

  for (chr in unique(as.character(runValue(seqnames(rowRanges(object))))))
  {
    genes <- genelist[seqnames(rowRanges(object)) == chr]
    grtrack <- Gviz::GeneRegionTrack(genes)
    atrack <- Gviz::AnnotationTrack(genes)
    gtrack <- Gviz::GenomeAxisTrack(add53 = TRUE)
    tracks <- c(gtrack, grtrack, atrack)
    z <- unlist(apply(scores(object), 2, function(u) {
      w <- unlist(as(u, "GRangesList"))
      lapply(split(w, strand(w), drop = TRUE), function(v)
        Gviz::DataTrack(v, type = "h", transformation = function(x) log2(x)))
    }
    ))
    
    Gviz::plotTracks(c(tracks, z), showId = TRUE)
  }
  
  par(oldpar)
}
