# test_seqcp-cumSeg.R
library(testthat)
context("seqcp algorithm callbacks")
library(GPosExperiment)
library(GenomicRanges)
library(rtracklayer)
rm(list = ls())

# test_that("create changepoints with cumSeg", {
#   set.seed(20180814)
#   devtools::load_all()
#   data("tdata_features")
#   data("tdata_gr_scores")
#   data("tdata_scores")
#
#   s1 <- GPosExperiment(features = tdata_features, scores = tdata_scores)
#   expect_equivalent(class(s1), "GPosExperiment", "new")
#   library(cumSeg)
#   parameters <- list(output = "3")
#   y <- calccp(s1,
#               algorithm = "cumSeg",
#               f = jumpoints,
#               parameters = parameters,
#               f2 = function(j) as.integer(j$psi))
#   result.list(s1) <- y
#   expect_equivalent(class(s1), "GPosExperiment", "after result.list<-")
#  save(s1, file = "../data/s1.rda")
#  z <- result.list(s1)
#  expect_equivalent(class(s1), "GPosExperiment", "after result.list<-+result.list")
#  expect_equivalent(class(z), "list")
#  expect_is(z[[1]], "seqcpSegments")
#
# TODO: DEBUG THIS
  # ax <- seqcp_from_files(feature_file = system.file("tests", "testfeatures_1.gff3",
  #                                                   package = "GPosExperiment", mustWork = TRUE),
  #                        scores_files = system.file("tests",
  #                                       paste0("testscores_1_", c("pos","neg"), ".bedGraph"),
  #                                                   package = "GPosExperiment", mustWork = TRUE),
  #                        sinfo = seqinfo(tdata_features),
  #                        id_column = "ID",
  #                        mcol_name = "score")
  # parameters <- list(output = "3")
  # ay <- calccp(ax,
  #             f = jumpoints,
  #             parameters = parameters,
  #             f2 = function(j) as.integer(j$psi))
  # result.list(ax) <- ay
  # expect_equivalent(result.list(x)[[1]]@results$changepoints, result.list(ax)[[1]]@results$changepoints)
# })
