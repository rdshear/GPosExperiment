features_file_1 <- system.file("extdata", "testfeatures_1.gff3.bgz", package = "GPosExperiment")
scores_file_1_pos <- system.file("extdata", "testscores_1_pos.bedGraph.bgz", package = "GPosExperiment")
scores_file_1_neg <- system.file("extdata", "testscores_1_neg.bedGraph.bgz", package = "GPosExperiment")

regions_of_interest <- import(features_file_1)
bedgraph_pos <- import(scores_file_1_pos)
strand(bedgraph_pos) <- "+"
bedgraph_neg <- import(scores_file_1_neg)
strand(bedgraph_neg) <- "-"
scores <- sort(c(bedgraph_pos, bedgraph_neg))
si <- seqinfo(regions_of_interest)


test_that("NETseqData simple constructor", {
  x <- NETseqData(sampleId = "S21", seqinfo = si)
  expect_s4_class(x, "NETseqData")
  expect_length(scores(x), 0)
  expect_length(subranges(x), 0)
  expect_equal(names(x), "S21")
})

test_that("NETseqData constructor with GRanges", {
  x <- NETseqData(sampleId = "x-15", scores = scores, subranges = regions_of_interest)
  expect_s4_class(x, "NETseqData")
  expect_length(scores(x), sum(width(scores)))
  expect_length(subranges(x), length(regions_of_interest))
  expect_equal(names(x), "x-15")
})

test_that("NETseqData constructor with stitched GPos", {
  w <- GPos(scores, stitch = TRUE)
  x <- NETseqData(sampleId = "x-15", scores = scores, subranges = regions_of_interest)
  expect_s4_class(x, "NETseqData")
  expect_length(scores(x), sum(width(scores)))
  expect_length(subranges(x), length(regions_of_interest))
  expect_equal(names(x), "x-15")
})

test_that("NEseqData constructor with explicit genome", {
  x <- NETseqData(seqinfo = Seqinfo(genome = "sacCer3"))
  expect_equivalent(genome(seqinfo(x))[1], "sacCer3")
  expect_equal(length(seqnames(seqinfo(x))), 17)
  expect_equal(seqnames(seqinfo(x))[17], "chrM")
})

test_that("NETseqDataFromBedgraph", {
  x <- NETseqDataFromBedgraph(sampleId = "xyz", 
        filename_pos = scores_file_1_pos, filename_neg = scores_file_1_neg, 
        seqinfo = si)
  expect_type(x, "list")
  expect_length(x, 1)
  x <- x[[1]]
  expect_s4_class(x, "NETseqData")
  total_reads <- sum(scores(x)$score)
  expected_reads <- sum(scores$score * width(scores))
  expect_equal(total_reads, expected_reads)
})

test_that("NETseqDataFromBedgraph two rows", {
  x <- NETseqDataFromBedgraph(sampleId = c("S1", "S2"),
                              filename_pos = rep(scores_file_1_pos, 2),
                              filename_neg = rep(scores_file_1_neg, 2),
                              filename_seg = rep(features_file_1, 2),
                              seqinfo = si)
  expect_type(x, "list")
  expect_length(x, 2)
  x <- x[[1]]
  expect_s4_class(x, "NETseqData")
  total_reads <- sum(scores(x)$score)
  expected_reads <- sum(scores$score * width(scores))
  expect_equal(total_reads, expected_reads)
  #TODO Better validity tests
})

