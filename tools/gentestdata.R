# gentestdata.R
library(BiocIO)
library(GenomicRanges)

set.seed(20180814)

# mini test genome
# TODOO: make a txdb
x_levels <- c("chrI", "chrII", "chrIII", "chrM")
x_lengths <-  c(1000, 2000, 3000, 500)
x_circ <- c(FALSE, FALSE, FALSE, TRUE)
x_genome <- "Xorg"

# MAKE FIRST TEST DATA
#features
tdata_features <- GRanges(c("chrI:20-250", "chrI:300-600", "chrII:1200-1600"))
names(tdata_features) <- c("ABC1", "XYZ2", "ABC")
seqlevels(tdata_features) <- x_levels
seqlengths(tdata_features) <- x_lengths
strand(tdata_features) <- c("+", "-", "+")
isCircular(tdata_features) <- x_circ
genome(tdata_features) <- x_genome
export(object = tdata_features,
       con = file.path("inst/extdata", "testfeatures_1.gff3"),
       version = "3",
       format = "GFF", index = TRUE)


step <- 73
tdata_scores <- lapply(c("+", "-"), function(w) RleList(lapply(c(chrI = 1000, chrII = 2000), function(u) {
  v <- as.vector(sapply(seq(from = 1, to = u + 1, by = step),
                        function(w) rep(w %% 11 + 1, step)))[1:u]
  Rle(rpois(v, lambda = v))
})))
names(tdata_scores) <- c("+", "-")
tdata_gr_scores <- GRangesList(lapply(tdata_scores, function(u) bindAsGRanges(score = u)))
for (i in names(tdata_gr_scores)) {
  strand(tdata_gr_scores[[i]]) <- i
}
for (i in names(tdata_gr_scores)) {
  f <- paste0("testscores_1_", ifelse(i == "+", "pos", "neg"), ".bedGraph")
  export(object = tdata_gr_scores[[i]],
         con = file.path("inst/extdata", f),
         format = "bedGraph", index = TRUE)
}

usethis::use_data(tdata_features, tdata_gr_scores, tdata_scores)
