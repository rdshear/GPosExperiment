test_that("NETseqData simple constructor", {
  refdata <- TestDataFilenames()
  sut <- NETseqData(sampleId = "S21", seqinfo = refdata$seqinfo)
  expect_s4_class(sut, "NETseqData")
  expect_length(scores(sut), 0)
  expect_length(subranges(sut), 0)
  expect_equal(names(sut), "S21")
})

test_that("NETseqData constructor with GRanges", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
  sut <- NETseqData(sampleId = sampleId, scores = test_bedgraphs$SRR12840066, 
                    subranges = genelist)
  expect_s4_class(sut, "NETseqData")
  expect_length(scores(sut), sum(width(test_bedgraphs$SRR12840066)))
  expect_length(subranges(sut), length(genelist))
  expect_equal(names(sut), sampleId)
})

test_that("NETseqData constructor with stitched GPos", {
  sampleId = "SRR12840066"
  refdata <- TestDataFilenames()
  w <- GPos(test_bedgraphs$SRR12840066, stitch = TRUE)
  sut <- NETseqData(sampleId = sampleId, scores = test_bedgraphs$SRR12840066, 
                    subranges = genelist)
  expect_s4_class(sut, "NETseqData")
  expect_length(scores(sut), sum(width(test_bedgraphs$SRR12840066)))
  expect_length(subranges(sut), length(genelist))
  expect_equal(names(sut), sampleId)
})

test_that("NEseqData constructor with explicit genome", {
  sut <- NETseqData(seqinfo = Seqinfo(genome = "sacCer3"))
  expect_equivalent(genome(seqinfo(sut))[1], "sacCer3")
  expect_equal(length(seqnames(seqinfo(sut))), 17)
  expect_equal(seqnames(seqinfo(sut))[17], "chrM")
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

test_that("NETseqDataFromBedgraph two rows", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
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
  # TODO: more verifcations
})
