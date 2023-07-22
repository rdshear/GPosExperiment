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

  # TODO split rownranges on seqname
  for (chr in unique(as.character(runValue(seqnames(rowRanges(object))))))
  {
    tracks <- list()
    genes <- rowRanges(object) # TODO WIP [seqnames(rowRanges(object)) == chr]
    tracks <- append(tracks, Gviz::GenomeAxisTrack(add53 = TRUE))
    tracks <- append(tracks, Gviz::GeneRegionTrack(genes,
                                                   name = chr,
                                                   shape = "arrow",
                                                   showId = TRUE,
                                                   symbol=genes$tx_name, 
                                                   add53 = TRUE))

    # TODO Add sample names to the DataTrack name: colNames(scores)
    z <- unlist(apply(scores(object), 2, function(u) {
      w <- unlist(as(u, "GRangesList"))
      w <- split(w, strand(w), drop=TRUE)
      mapply(function(v, id)
        Gviz::DataTrack(v, type = "h", name=id,
                        transformation = function(x) log2(x)),
        v = w, id = names(w))
    }))
        
    
    Gviz::plotTracks(c(tracks, z))
  }
  
  par(oldpar)
}
