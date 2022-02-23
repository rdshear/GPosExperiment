test_that("Instantion of GPosExperiement 11rx1c", {
  refdata <- TestDataFilenames()
  sampleId = "SRR12840066"
  nsd <-
    NETseqData(scores = test_bedgraphs$SRR12840066, sampleId = sampleId)
  sut <- GPosExperiment(sample = nsd, rowRanges = genelist,
                        seqinfo = refdata$seqinfo)
  s <- scores(sut)
  expect_true(all(width(genelist) == sapply(s, length)))
})

test_that("Instantion of GPosExperiement 11rx2c", {
  refdata <- TestDataFilenames()
  nsd <- mapply(function(s, id) {
    NETseqData(scores = s, sampleId = id,
               seqinfo = refdata$seqinfo)
  }, s = test_bedgraphs, id = names(test_bedgraphs))
  sut <- GPosExperiment(sample = nsd,
                        rowRanges = genelist)
  expect_s4_class(sut, "GPosExperiment")
  s <- scores(sut)
  sl <- sapply(s, length)
  dim(sl) <- dim(s)

  expect_true(all(sl[,1] == sl[,2]))
  expect_true(all(width(genelist) == sl[, 1]))
  # TODO more expect_*
  })
