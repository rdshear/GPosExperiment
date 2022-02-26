test_that("Instantion of GPosExperiement 11rx1c", {
  refdata <- .TestDataFilenames()
  sampleId = "SRR12840066"
  nsd <-
    NETseqData(scores = test_bedgraphs$SRR12840066, 
               mask = test_masks$SRR12840066,
               sampleId = sampleId)
  sut <- GPosExperiment(sample = nsd, rowRanges = genelist,
                        seqinfo = refdata$seqinfo)
  s <- vscores(sut)
  expect_true(all(width(genelist) == sapply(s, length)))
  m <- mask(sut)
  mask_widths <- sapply(m, function(u) sum(width(u)))
  dim(mask_widths) <- dim(m)
  dimnames(mask_widths) <- dimnames(m)
  expect_equal(sum(mask_widths), sum(width(test_masks$SRR12840066)))
})

test_that("Instantion of GPosExperiement 11rx2c", {
  refdata <- .TestDataFilenames()
  nsd <- mapply(function(s, m, id) {
    NETseqData(scores = s, sampleId = id, mask = m,
               seqinfo = refdata$seqinfo)
  }, s = test_bedgraphs, m = test_masks, id = names(test_bedgraphs))
  sut <- GPosExperiment(sample = nsd,
                        rowRanges = genelist)
  expect_s4_class(sut, "GPosExperiment")
  s <- vscores(sut)
  sl <- sapply(s, length)
  dim(sl) <- dim(s)

  expect_true(all(sl[,1] == sl[,2]))
  expect_true(all(width(genelist) == sl[, 1]))
  # TODO more expect_*
  })
