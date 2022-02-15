# test_GPosExpSE-cumSeg.R
library(testthat)
context("GPosExperiment plot")
library(GPosExperiment)
library(GenomicRanges)
library(rtracklayer)
# rm(list = ls())
# 
# # TODO: debug result.list(s1) <- y problem
# test_that("plot changepoints", {
#   set.seed(20180814)
#   devtools::load_all()
#   data("tdata_features")
#   data("tdata_gr_scores")
#   data("tdata_scores")
# 
#   s1 <- GPosExperiment(features = tdata_features, scores = tdata_scores)
#   expect_equivalent(class(s1), "GPosExperiment", "new")
# 
#   ##### cumSeg
#   library(cumSeg, quietly = TRUE)
#   parameters <- list(output = "3")
#   y <- calccp(s1, "cumSeg",
#               f = jumpoints,
#               parameters = parameters,
#               f2 = function(j) as.integer(j$psi))
#   result.list(s1) <- y
#   expect_equivalent(class(s1), "GPosExperiment", "after cumSeg")
# 
#   ##### Segmentor3IsBack
#   library(Segmentor3IsBack, quietly = TRUE)
# 
#   # TODO: move the mcparallel to package...invoke with new parameter isolation.mode = TRUE
#   get_segmentation <- function(score) {
#     p <- parallel:::mcparallel(get_seg_inner(score))
#     r <- parallel:::mccollect(p)
#     r[[1]]
#   }
# 
#   # HACK: Hard coded Kmax
#   get_seg_inner <- function(score) {
#     seg <- Segmentor(score, Kmax = 6, model = 3, keep = TRUE)
#     kChoose <- SelectModel(seg, penalty = "BIC")
#     a <- getBreaks(seg)[kChoose, ]
#     # trim zeros
#     a[a > 0]
#   }
# 
#   # TODO: refactor get_segmentation & get_seg_inner to use parameters
#   parameters <- list()
# 
#   y <- calccp(s1, "Segmentor3IsBack",
#               f = get_segmentation,
#               parameters = parameters,
#               f2 = function(j) j)
#   result.list(s1) <- y
#   expect_equivalent(class(s1), "GPosExperiment", "after Segmentor3IsBack")
# 
#   ##### breakpoint
#   library(breakpoint)
# 
#   parameters <- list(Nmax = s1@kMax, parallel = TRUE)
#   # TODO: Skipping this due to long run times
#   # y <- calccp(s1, "breakpoint",
#   #             f = function(score, ...) CE.NB(data.frame(score), ...),
#   #             parameters = parameters,
#   #             f2 = function(seg) seg$BP.Loc)
#   # result.list(s1) <- y
#   # expect_equivalent(class(s1), "GPosExperiment", "after breakpoint")
# 
#   ############   EBS
#     library(EBS)
#   # TODO: Skipping this due to long run times
#     #
#     # y <- calccp(s1, "EBS",
#     #        f = EBSegmentation,
#     #        parameters <- list(model = 3, Kmax = s1@kMax),
#     #        f2 = function(seg) {
#     #          icl <- EBSICL(seg)
#     #          k <- icl$NbICL
#     #          breaks <-
#     #            sapply(seq(k - 1), function(u)
#     #              which.max(EBSDistrib(seg, u, k)))
#     #          breaks
#     #        })
#     # # TODO: replace results.list replacement function with addtional logic in calccp
#     # result.list(s1) <- y
#     # expect_equivalent(class(s1), "GPosExperiment", "after EBS")
# 
#     # TODO: make summary method for GPosExperiment that gives counts of errors, etc
#   #plot(s1, filter = names(s1@features), tracks = seq_along(s1@result.list))
#     plot(s1)
# })
