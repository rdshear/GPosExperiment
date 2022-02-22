test_that("Instantion of GPosExperiement Nrx1c", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
  nsd <-
    NETseqData(scores = test_bedgraphs$SRR12840066, sampleId = sampleId)
  sut <- GPosExperiment(sample = nsd, rowRanges = genelist,
                        seqinfo = refdata$seqinfo)
  s <- scores(sut)
  # TODO: need to deal with sparse matrix
  expect_true(all(width(genelist) == sapply(s, length)))
})

test_that("Instantion of GPosExperiement Nrx3c", {
  refdata <- TestDataFilenames()
  nsd <- mapply(function(s, id) {
    NETseqData(scores = s, sampleId = id,
               seqinfo = refdata$seqinfo)
  }, s = test_bedgraphs, id = names(test_bedgraphs))
  sut <- GPosExperiment(sample = nsd,
                        rowRanges = genelist)
  expect_s4_class(sut, "GPosExperiment")
  # TODO more expect_*
  })
