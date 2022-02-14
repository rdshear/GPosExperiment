library(testthat)
context("GPosExperiment basics")

library(rtracklayer)
library(SummarizedExperiment)
library(GPosExperiment)

rm(list = ls())
set.seed(20180814)

#' TODO: Set up test data: in inst/extdata directory
#' reference per following:  system.file("extdata", "SampleTestFile.txt", package = "GPosExperiment")

features_file_1 <- system.file("extdata", "testfeatures_1.gff3.bgz", package = "GPosExperiment")
scores_file_1_pos <- system.file("extdata", "testscores_1_pos.bedGraph.bgz", package = "GPosExperiment")
scores_file_1_neg <- system.file("extdata", "testscores_1_neg.bedGraph.bgz", package = "GPosExperiment")

test_that("Basic instantiation of GPosExperiment object", {
  x <- GPosExperiment()
  expect_s4_class(x, "GPosExperiment")
})

test_that("Basic instantiation of GPosExpSamples object", {
  x <- GPosExpSamples()
  expect_s4_class(x, "GPosExpSamples")
})

test_that("GPosExpSamples object with filename references", {
  x <- GPosExpSamples(scoreFileDirectory = system.file("extdata", package = "GPosExperiment"),
                   segmentFileDirectory  = system.file("extdata", package = "GPosExperiment"),
                   sampleNames = "S01",
                   scoreFilesPos = "testscores_1_pos.bedGraph.bgz",
                   scoreFilesNeg = "testscores_1_neg.bedGraph.bgz",
                   segmentFiles = "testfeatures_1.gff3.bgz")

  expect_equivalent(class(x), "GPosExpSamples")
})

test_that("GPosExpSamples object with unequal sample vector lengths", {
  expect_error({
    GPosExpSamples(scoreFileDirectory = system.file("extdata", package = "GPosExperiment"),
                segmentFileDirectory  = system.file("extdata", package = "GPosExperiment"),
                sampleNames = c("S01", "T001"),
                scoreFilesPos = "xtestscores_1_pos.bedGraph.bgz",
                scoreFilesNeg = "testscores_1_neg.bedGraph.bgz",
                segmentFiles = "testfeatures_1.gff3.bgz")
  }, message = "invalid class “GPosExpSamples” object: (+) strand score file not found")
})
test_that("GPosExpSamples object with missing file", {
  expect_error({
    GPosExpSamples(scoreFileDirectory = system.file("extdata", package = "GPosExperiment"),
                 segmentFileDirectory  = system.file("extdata", package = "GPosExperiment"),
                 sampleNames = "S01",
                 scoreFilesPos = "xtestscores_1_pos.bedGraph.bgz",
                 scoreFilesNeg = "testscores_1_neg.bedGraph.bgz",
                 segmentFiles = "testfeatures_1.gff3.bgz")
  }, message = "invalid class “GPosExpSamples” object: (+) strand score file not found")
})
#
# test_that("Create object with scores from GRanges object", {
#   x <- GPosExperiment(tdata_features, tdata_scores)
#   expect_equivalent(class(x), "GPosExperiment")
# })

# TODO: Tests for all GPosExperimentSegments slots
# TODO: Tests for all GPosExperiment slots
# TODO: reinstate when scores <- is fixed
# test_that("Set scores", {
#   x <- GPosExperiment(tdata_features, scores = tdata_scores)
#   expect_equivalent(class(x), "GPosExperiment")
#   y <- GPosExperiment(tdata_features)
#   scores(y) <- tdata_gr_scores
#   expect_identical(x, y)
# })

# TODO: Test inverse function for GRanges vs RleList
# TODO: Test trim /fill out of scores vs seqinfo

# test_that("Create GPosExperimentSegments", {
#   x <- GPosExperiment(tdata_features, scores = tdata_scores)
#   x <- GPosExperimentSegments("cumSeq")
#   expect_equivalent(class(x), "GPosExperimentSegments")
#   expect_identical(x@parameters, list())
#   expect_identical(x@results, data.frame())
# })