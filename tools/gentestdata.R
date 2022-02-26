# gentestdata.R
library(rtracklayer)
library(GenomicRanges)

set.seed(20180814)

## data from datasets
# TODO Document your data (see 'https://r-pkgs.org/data.html')
a <- .TestDataFilenames()
genelist <- import(a$genes)
names(genelist) <- genelist$ID
seqinfo(genelist) <- a$seqinfo
usethis::use_data(genelist)

test_bedgraphs <- apply(a$samples, 1, function(x){
  pos <- import(x["bedgraph_pos"])
  strand(pos) <- "+"
  neg <- import(x["bedgraph_neg"])
  strand(neg) = "-"
  y <- GenomicRanges::sort(c(pos, neg))
  seqinfo(y) <- a$seqinfo
  y
})

use_data(test_bedgraphs)

test_masks <- apply(a$samples, 1, function(u) {
  m <- import(u["mask"])
  seqinfo(m) <- a$seqinfo
  m
})

use_data(test_masks)
