test_that("Instantion of GPosExperiement Nrx1c", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
  nsd <-
    NETseqData(scores = test_bedgraphs$SRR12840066, sampleId = sampleId)
  sut <- GPosExperiment(sample = nsd, rowRanges = genelist,
                        seqinfo = refdata$seqinfo)
  s <- scores(sut)
  # TODO: need to deal with sparce matrix
  expect_true(all(width(genelist) == sapply(s, length)))
})

test_that("Instantion of GPosExperiement 2rx3c", {
  refdata <- TestDataFilenames()
  nsd <- mapply(function(s, id) {
    NETseqData(scores = s, sampleId = id,
               seqinfo = refdata$seqinfo)
  }, s = test_bedgraphs, id = names(test_bedgraphs))
  sut <- GPosExperiment(sample = nsd,
                        rowRanges = genelist)
  expect_s4_class(sut, "GPosExperiment")
  })


test_that("Double instantiation of GPosExperiment from Bedgraph files", {
  nsd <- NETseqDataFromBedgraph(
    sampleId = c("S01", "T01"),
    filename_pos = rep(scores_file_1_pos, 2),
    filename_neg = rep(scores_file_1_neg, 2),
    filename_seg = rep(features_file_1, 2),
    seqinfo = si
  )
  sut <-
    GPosExperiment(rowRanges = ff_ranges,
                   sample = nsd,
                   seqinfo = si)
  
  expect_s4_class(sut, "GPosExperiment")
})
