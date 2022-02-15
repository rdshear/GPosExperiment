library(GenomicRanges)
library(rtracklayer)

features_file_1 <- system.file("extdata", "testfeatures_1.gff3.bgz", package = "GPosExperiment")
scores_file_1_pos <- system.file("extdata", "testscores_1_pos.bedGraph.bgz", package = "GPosExperiment")
scores_file_1_neg <- system.file("extdata", "testscores_1_neg.bedGraph.bgz", package = "GPosExperiment")

regions_of_interest <- import(features_file_1)
bedgraph_pos <- import(scores_file_1_pos)
strand(bedgraph_pos) <- "+"
bedgraph_neg <- import(scores_file_1_neg)
strand(bedgraph_neg) <- "-"
scores <- sort(c(bedgraph_pos, bedgraph_neg))


test_that("NETseqData simple constructor", {
  x <- NETseqData(sampleId = "S21")
  expect_s4_class(x, "NETseqData")
  expect_length(x@scores, 0)
  expect_length(x@segments, 0)
  expect_equal(x@sampleId, "S21")
})

test_that("NETseqData constructor with GRanges", {
  x <- NETseqData(sampleId = "x-15", scores = scores, segments = regions_of_interest)
  expect_s4_class(x, "NETseqData")
  expect_length(x@scores, sum(width(scores)))
  expect_length(x@segments, length(regions_of_interest))
  expect_equal(x@sampleId, "x-15")
})

test_that("NETseqData constructor with stitched GPos", {
  w <- GPos(scores, stitch = TRUE)
  x <- NETseqData(sampleId = "x-15", scores = scores, segments = regions_of_interest)
  expect_s4_class(x, "NETseqData")
  expect_length(x@scores, sum(width(scores)))
  expect_length(x@segments, length(regions_of_interest))
  expect_equal(x@sampleId, "x-15")
})
