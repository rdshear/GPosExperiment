library(rtracklayer)
features_file_1 <- system.file("extdata", "testfeatures_1.gff3.bgz", package = "GPosExperiment")
scores_file_1_pos <- system.file("extdata", "testscores_1_pos.bedGraph.bgz", package = "GPosExperiment")
scores_file_1_neg <- system.file("extdata", "testscores_1_neg.bedGraph.bgz", package = "GPosExperiment")
si <- GenomeInfoDb::Seqinfo(genome = "sacCer3")
ff_ranges <- import(features_file_1)
names(ff_ranges) <- ff_ranges$ID
seqinfo(ff_ranges) <- si



test_that("Basic instantiation of GPosExperiment object", {
  nsdf <- NETseqDataFromBedgraph(sampleId = c("S01", "T01"),  
                               filename_pos = rep(scores_file_1_pos, 2),
                               filename_neg = rep(scores_file_1_neg, 2),
                               filename_seg = rep(features_file_1, 2),
                               seqinfo = si)
  x <- GPosExperiment(rowRanges = ff_ranges, sample = nsdf, seqinfo = si)
  expect_s4_class(x, "GPosExperiment")
})

