test_that("HRseqData simple constructor", {
  refdata <- .TestDataFilenames()
  # TODO add at least one gene from a different chromosome
  # TODO fix up samples default
  sut <- HRseqData(sampleId = "S21", 
                    scores = GRanges(seqinfo = refdata$seqinfo), 
                    seqinfo = refdata$seqinfo)
  expect_s4_class(sut, "HRseqData")
  # TODO function return inconsistency, scores here is a GPos, and for the SE it is 
  # integer list
  expect_length(scores(sut), 0)
  expect_length(mask(sut), 0)
  expect_equal(names(sut), "S21")
})

test_that("HRseqData constructor with GRanges", {
  refdata <- .TestDataFilenames()
  sampleId = "SRR12840066"
  sut <- HRseqData(sampleId = sampleId, scores = test_bedgraphs$SRR12840066, 
                    mask = genelist[2,3])
  expect_s4_class(sut, "HRseqData")
  expect_equal(sum(scores(sut)$score), 
               sum(width(test_bedgraphs$SRR12840066) * test_bedgraphs$SRR12840066$score))
  expect_equal(mask(sut), genelist[2,3])
  expect_equal(names(sut), sampleId)
})

test_that("HRseqData constructor with stitched GPos", {
  sampleId = "SRR12840066"
  refdata <- .TestDataFilenames()
  w <- GPos(test_bedgraphs$SRR12840066, stitch = TRUE)
  sut <- HRseqData(sampleId = sampleId, scores = test_bedgraphs$SRR12840066, 
                    mask = genelist[2])
  expect_s4_class(sut, "HRseqData")
  expect_equal(sum(scores(sut)$score), 
               sum(width(test_bedgraphs$SRR12840066) * 
                     test_bedgraphs$SRR12840066$score))
  expect_equal(names(sut), sampleId)
})

test_that("HRseqDataFromBedgraph", {
  refdata <- .TestDataFilenames()
  sampleId = "SRR12840066"
  sampleParameters <- refdata$samples[sampleId, ]
  sut <- HRseqDataFromBedgraph(sampleId = sampleId, 
        filename_pos = sampleParameters$bedgraph_pos, 
        filename_neg = sampleParameters$bedgraph_neg,
        seqinfo = refdata$seqinfo)
  expect_type(sut, "list")
  expect_length(sut, 1)
  sut <- sut[[1]]
  expect_s4_class(sut, "HRseqData")
  scores <- scores(sut)
  total_reads <- sum(scores$score)
  ref <- test_bedgraphs$SRR12840066
  expected_reads <- sum(ref$score * width(ref))
  expect_equal(total_reads, expected_reads)
  zero_filled_scores <- scores(sut, zero_fill = TRUE)
  expect_equal(total_reads, sum(zero_filled_scores$score))
  # we know that our test run has only scores in chrI:+ and chrI:-
  genome_size <- sum(seqlengths(seqinfo(sut)["chrI"]))
  expect_equal(genome_size * 2, length(zero_filled_scores))
  r <- GRanges("chrI:100-500:+")
  ref_scoped <- .OverlappedRanges(r, ref)
  scores_scoped <- scores(sut, range_scope = r)
  expect_equal(sum(scores_scoped$score), 
               sum(ref_scoped$score * width(ref_scoped)))
  zscores_scoped <- scores(sut, range_scope = r, zero_fill = TRUE)
  expect_equal(sum(zscores_scoped$score), 
               sum(ref_scoped$score * width(ref_scoped)))
})

test_that("HRseqDataFromBedgraph two samples", {
  sampleId = "SRR12840066"
  refdata <- .TestDataFilenames()
  sampledata <- refdata$samples
  sut <- HRseqDataFromBedgraph(sampleId = sampledata$sampleId,
                              filename_pos = sampledata$bedgraph_pos,
                              filename_neg = sampledata$bedgraph_neg,
                              seqinfo = refdata$seqinfo)
  expect_type(sut, "list")
  expect_length(sut, 2)
  for (u in sut) {
    expect_s4_class(u, "HRseqData")
  }
  total_reads <- sum(scores(sut[[sampleId]])$score)
  expected_reads <- sum(test_bedgraphs$SRR12840066$score * width(test_bedgraphs$SRR12840066))
  expect_equal(total_reads, expected_reads)
  
})

test_that("HRseqData from GRanges scores", {
  nsd <-HRseqData(scores = test_bedgraphs$SRR12840066, sampleId = "SRR12840066")
  reads_from_granges <- sum(width(test_bedgraphs$SRR12840066) * test_bedgraphs$SRR12840066$score)
  reads_from_HRseqData <- sum(scores(nsd)$score)
  expect_equal(reads_from_HRseqData, reads_from_granges)
})


test_that("HRseqData from BAM file",{
  refdata <- .TestDataFilenames()

  sut <- HRseqDataFromBAM(refdata$samples$sampleId,
                           refdata$samples$bam,
                           genelist,
                           refdata$seqinfo)
  
  expect_type(sut, "list")
  expect_true(all(sapply(sut, isa,"HRseqData")))
  for (s in seq_along(sut)) {
    expect_equal(as.character(mask(sut[[s]])), 
               as.character(test_masks[[s]]))
  }
  # TODO Break into 2 tests...after setting up appropriate fixture
  
  nsd <- sut
  sut <- GPosExperiment(rowRanges = genelist, sample = nsd, seqinfo = si)
  # 
  expect_s4_class(sut, "GPosExperiment")
  s <- scores(sut)
  sc <- apply(s, 1:2, function(u) sum(u[[1]]$score, na.rm = TRUE))
  sna <- apply(s, 1:2, function(u) sum(is.na(u[[1]]$score)))
  
  snm <- scores(sut, apply_mask = FALSE)
  snm_c <- apply(snm, 1:2, function(u) sum(u[[1]]$score, na.rm = TRUE))
  snm_na <- apply(snm, 1:2, function(u) sum(is.na(u[[1]]$score)))
  expect_true(all(sc <= snm_c))
  # TODO: match sna vs number of positions lost to dups
  # TODO: zero_fill = TRUE
  # TODO: more verification
})
