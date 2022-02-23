test_that("NETseqData simple constructor", {
  refdata <- TestDataFilenames()
  # TODO add at least one gene from a different chromosome
  # TODO fix up samples default
  sut <- NETseqData(sampleId = "S21", 
                    scores = GRanges(seqinfo = refdata$seqinfo), 
                    seqinfo = refdata$seqinfo)
  expect_s4_class(sut, "NETseqData")
  # TODO function return inconsistency, scores here is a GPos, and for the SE it is 
  # integer list
  expect_length(scores(sut), 0)
  expect_length(mask(sut), 0)
  expect_equal(names(sut), "S21")
})

test_that("NETseqData constructor with GRanges", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
  sut <- NETseqData(sampleId = sampleId, scores = test_bedgraphs$SRR12840066, 
                    mask = genelist[2,3])
  expect_s4_class(sut, "NETseqData")
  expect_equal(sum(scores(sut)$score), 
               sum(width(test_bedgraphs$SRR12840066) * test_bedgraphs$SRR12840066$score))
  expect_equal(mask(sut), genelist[2,3])
  expect_equal(names(sut), sampleId)
})

test_that("NETseqData constructor with stitched GPos", {
  sampleId = "SRR12840066"
  refdata <- TestDataFilenames()
  w <- GPos(test_bedgraphs$SRR12840066, stitch = TRUE)
  sut <- NETseqData(sampleId = sampleId, scores = test_bedgraphs$SRR12840066, 
                    mask = genelist[2])
  expect_s4_class(sut, "NETseqData")
  expect_equal(sum(scores(sut)$score), 
               sum(width(test_bedgraphs$SRR12840066) * test_bedgraphs$SRR12840066$score))
  expect_equal(names(sut), sampleId)
})

test_that("NETseqDataFromBedgraph", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
  sampleParameters <- refdata$samples[sampleId, ]
  sut <- NETseqDataFromBedgraph(sampleId = sampleId, 
        filename_pos = sampleParameters$bedgraph_pos, 
        filename_neg = sampleParameters$bedgraph_neg,
        seqinfo = refdata$seqinfo)
  expect_type(sut, "list")
  expect_length(sut, 1)
  sut <- sut[[1]]
  expect_s4_class(sut, "NETseqData")
  total_reads <- sum(scores(sut)$score)
  expected_reads <- sum(test_bedgraphs$SRR12840066$score * width(test_bedgraphs$SRR12840066))
  expect_equal(total_reads, expected_reads)
})

test_that("NETseqDataFromBedgraph two samples", {
  sampleId = "SRR12840066"
  refdata <- TestDataFilenames()
  sampledata <- refdata$samples
  sut <- NETseqDataFromBedgraph(sampleId = sampledata$sampleId,
                              filename_pos = sampledata$bedgraph_pos,
                              filename_neg = sampledata$bedgraph_neg,
                              seqinfo = refdata$seqinfo)
  expect_type(sut, "list")
  expect_length(sut, 2)
  for (u in sut) {
    expect_s4_class(u, "NETseqData")
  }
  total_reads <- sum(scores(sut[[sampleId]])$score)
  expected_reads <- sum(test_bedgraphs$SRR12840066$score * width(test_bedgraphs$SRR12840066))
  expect_equal(total_reads, expected_reads)
  
})

test_that("NETseqData from GRanges scores", {
  nsd <-NETseqData(scores = test_bedgraphs$SRR12840066, sampleId = "SRR12840066")
  reads_from_granges <- sum(width(test_bedgraphs$SRR12840066) * test_bedgraphs$SRR12840066$score)
  reads_from_NETseqData <- sum(scores(nsd)$score)
  expect_equal(reads_from_NETseqData, reads_from_granges)
})


test_that("NETseqData from BAM file",{
  refdata <- TestDataFilenames()

  sut <- NETseqDataFromBAM(refdata$samples$sampleId,
                           refdata$samples$bam,
                           genelist,
                           refdata$seqinfo)
  
  expect_type(sut, "list")
  expect_true(all(sapply(sut, isa,"NETseqData")))
  # TODO Break into 2 tests...after setting up appropriate fixture
  
  nsd <- sut
  sut <- GPosExperiment(rowRanges = genelist, sample = nsd, seqinfo = si)
  # 
  expect_s4_class(sut, "GPosExperiment")
  s <- vscores(sut)
  sc <- sapply(s, sum, na.rm = TRUE)
  dim(sc) <- dim(s)
  sna <- sapply(s, function(u) sum(is.na(u)))
  dim(sna) <- dim(s)
  
  snm <- vscores(sut, apply.mask = FALSE)
  snm_c <- sapply(snm, sum, na.rm = TRUE)
  dim(snm_c) <- dim(snm)
  snm_na <- sapply(snm, function(u) sum(is.na(u)))
  dim(snm_na) <- dim(snm)
  expect_true(all(sc <= snm_c))
  # TODO: match sna vs number of positions lost to dups
  # TODO: more verification
})
