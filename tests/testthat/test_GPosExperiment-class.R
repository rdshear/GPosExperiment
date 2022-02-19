library(rtracklayer)
features_file_1 <- system.file("extdata", "testfeatures_1.gff3.bgz", package = "GPosExperiment")
scores_file_1_pos <- system.file("extdata", "testscores_1_pos.bedGraph.bgz", package = "GPosExperiment")
scores_file_1_neg <- system.file("extdata", "testscores_1_neg.bedGraph.bgz", package = "GPosExperiment")
si <- GenomeInfoDb::Seqinfo(genome = "sacCer3")
ff_ranges <- import(features_file_1)
names(ff_ranges) <- ff_ranges$ID
seqinfo(ff_ranges) <- si
data("tdata_features")
data("tdata_scores")
data("tdata_gr_scores")

test_that("Instantion of GPosExperiement 3rx1c", {
  nsd <-NETseqData(scores = unlist(tdata_gr_scores), sampleId = "S01")
  sut <- GPosExperiment(sample = nsd, rowRanges = tdata_features)
  s <- scores(sut)
  # TODO: need to deal with sparce matrix
  expect_true(all(width(tdata_features) == sapply(s, length)))
})

test_that("Instantion of GPosExperiement 3rx3c", {
  nsd <-NETseqData(scores = unlist(tdata_gr_scores), sampleId = "S01")
  nsd1 <- NETseqData(scores = tdata_gr_scores[[1]], sampleId = "S02")
  nsd2 <- NETseqData(scores = tdata_gr_scores[[2]], sampleId = "S03")
  sut <- GPosExperiment(sample = list(nsd, nsd1, nsd2), 
                      rowRanges = tdata_features)
  s <- scores(sut)
  v <- sapply(s, length)
  dim(v) <- dim(s)
  u <- matrix(v[,1], ncol = 3, nrow = 3)
  u[2,2] <- 0
  u[1,3] <- 0
  u[3,3] <- 0
  expect_equal(u, v)
})


test_that("Double instantiation of GPosExperiment from Bedgraph files", {
  nsd <- NETseqDataFromBedgraph(sampleId = c("S01", "T01"),  
                               filename_pos = rep(scores_file_1_pos, 2),
                               filename_neg = rep(scores_file_1_neg, 2),
                               filename_seg = rep(features_file_1, 2),
                               seqinfo = si)
  sut <- GPosExperiment(rowRanges = ff_ranges, sample = nsd, seqinfo = si)

  expect_s4_class(sut, "GPosExperiment")
})

