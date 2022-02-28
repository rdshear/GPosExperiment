test_that("Instantion of GPosExperiement 11rx1c", {
  refdata <- .TestDataFilenames()
  sampleId = "SRR12840066"
  nsd <-
    NETseqData(scores = test_bedgraphs$SRR12840066, 
               mask = test_masks$SRR12840066,
               sampleId = sampleId)
  sut <- GPosExperiment(sample = nsd, rowRanges = genelist,
                        seqinfo = refdata$seqinfo)
  s <- scores(sut)
  # expect_true(all(width(genelist) == sapply(s, length)))
  m <- mask(sut)
  mask_widths <- .matrix_apply(m, function(u) sum(width(u)))
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
  s <- scores(sut)
  sl <- .matrix_apply(s, length)

  expect_true(all(sl[,1] == sl[,2]))
  expect_true(all(width(genelist) == sl[, 1]))
  # TODO more expect_*
  })
