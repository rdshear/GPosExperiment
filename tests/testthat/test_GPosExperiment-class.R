test_that("Instantion of GPosExperiement 11rx1c", {
  refdata <- .TestDataFilenames()
  sampleId = "SRR12840066"
  refmasks <- test_masks$SRR12840066
  nsd <-
    HRseqData(scores = test_bedgraphs$SRR12840066, 
               mask = refmasks,
               sampleId = sampleId)
  sut <- HRseqExperiment(sample = nsd, rowRanges = genelist,
                        seqinfo = refdata$seqinfo)

  ref_scores <- .OverlappedRanges(genelist, test_bedgraphs$SRR12840066)
  ref_totalreads <- sum(ref_scores$score * width(ref_scores))
  s <- scores(sut, zero_fill = FALSE, apply_mask = FALSE)
  sl <- apply(s, 1:2, function(u) length(u[[1]]))
  sc <- apply(s, 1:2, function(u) sum(u[[1]]$score))
  expect_true(all(sl <= width(genelist)))
  expect_equal(ref_totalreads, sum(sc))
  
  sz <- scores(sut, zero_fill = TRUE, apply_mask = FALSE)
  szl <- apply(sz, 1:2, function(u) length(u[[1]]))
  szc <- apply(sz, 1:2, function(u) sum(u[[1]]$score))
  expect_true(all(szl == width(genelist)))
  expect_equal(ref_totalreads, sum(szc))
  # expect_true(all(width(genelist) == sapply(s, length)))
  m <- mask(sut)
  mask_widths <- apply(m, 1:2, function(u) sum(width(u[[1]])))
  expect_equal(sum(mask_widths), sum(width(refmasks)))
  # scores <- scores(sut)
  # total_reads <- sum(scores$score)
  # ref <- test_bedgraphs$SRR12840066
  # expected_reads <- sum(ref$score * width(ref))
  # expect_equal(total_reads, expected_reads)
  # zero_filled_scores <- scores(sut, zero_fill = TRUE)
  # expect_equal(total_reads, sum(zero_filled_scores$score))
  # genome_size <- sum(seqlengths(seqinfo(sut)))
  # expect_equal(genome_size * 2, length(zero_filled_scores))
  # r <- GRanges("chrI:100-500:+")
  # ref_scoped <- .OverlappedRanges(r, ref)
  # scores_scoped <- scores(sut, range_scope = r)
  # expect_equal(sum(sum(scores_scoped$score)), 
  #              sum(ref_scoped$score * width(ref_scoped)))
})

test_that("Instantion of GPosExperiement 11rx2c", {
  refdata <- .TestDataFilenames()
  nsd <- mapply(function(s, m, id) {
    HRseqData(scores = s, sampleId = id, mask = m,
               seqinfo = refdata$seqinfo)
  }, s = test_bedgraphs, m = test_masks, id = names(test_bedgraphs))
  sut <- HRseqExperiment(sample = nsd,
                        rowRanges = genelist)
  expect_s4_class(sut, "HRseqExperiment")
  s <- scores(sut, zero_fill = TRUE, apply_mask = FALSE)
  sl <- apply(s, 1:2, function(u) length(u[[1]]))

  expect_true(all(sl[,1] == sl[,2]))
  expect_true(all(width(genelist) == sl[, 1]))
  # TODO more expect_*
  })
