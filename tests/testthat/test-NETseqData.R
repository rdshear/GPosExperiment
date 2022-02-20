features_file_1 <- system.file("extdata", "testfeatures_1.gff3.bgz", package = "GPosExperiment")
scores_file_1_pos <- system.file("extdata", "testscores_1_pos.bedGraph.bgz", package = "GPosExperiment")
scores_file_1_neg <- system.file("extdata", "testscores_1_neg.bedGraph.bgz", package = "GPosExperiment")
data("tdata_features")
data("tdata_scores")
data("tdata_gr_scores")

regions_of_interest <- import(features_file_1)
bedgraph_pos <- import(scores_file_1_pos)
strand(bedgraph_pos) <- "+"
bedgraph_neg <- import(scores_file_1_neg)
strand(bedgraph_neg) <- "-"
scores <- sort(c(bedgraph_pos, bedgraph_neg))
si <- seqinfo(regions_of_interest)


test_that("NETseqData simple constructor", {
  sut <- NETseqData(sampleId = "S21", seqinfo = si)
  expect_s4_class(sut, "NETseqData")
  expect_length(scores(sut), 0)
  expect_length(subranges(sut), 0)
  expect_equal(names(sut), "S21")
})

test_that("NETseqData constructor with GRanges", {
  sut <- NETseqData(sampleId = "X-15", scores = scores, subranges = regions_of_interest)
  expect_s4_class(sut, "NETseqData")
  expect_length(scores(sut), sum(width(scores)))
  expect_length(subranges(sut), length(regions_of_interest))
  expect_equal(names(sut), "X-15")
})

test_that("NETseqData constructor with stitched GPos", {
  w <- GPos(scores, stitch = TRUE)
  sut <- NETseqData(sampleId = "X-15", scores = scores, subranges = regions_of_interest)
  expect_s4_class(sut, "NETseqData")
  expect_length(scores(sut), sum(width(scores)))
  expect_length(subranges(sut), length(regions_of_interest))
  expect_equal(names(sut), "X-15")
})

test_that("NEseqData constructor with explicit genome", {
  sut <- NETseqData(seqinfo = Seqinfo(genome = "sacCer3"))
  expect_equivalent(genome(seqinfo(sut))[1], "sacCer3")
  expect_equal(length(seqnames(seqinfo(sut))), 17)
  expect_equal(seqnames(seqinfo(sut))[17], "chrM")
})

test_that("NETseqDataFromBedgraph", {
  sut <- NETseqDataFromBedgraph(sampleId = "xyz", 
        filename_pos = scores_file_1_pos, filename_neg = scores_file_1_neg, 
        seqinfo = si)
  expect_type(sut, "list")
  expect_length(sut, 1)
  sut <- sut[[1]]
  expect_s4_class(sut, "NETseqData")
  total_reads <- sum(scores(sut)$score)
  expected_reads <- sum(scores$score * width(scores))
  expect_equal(total_reads, expected_reads)
})

test_that("NETseqDataFromBedgraph two rows", {
  sut <- NETseqDataFromBedgraph(sampleId = c("S1", "S2"),
                              filename_pos = rep(scores_file_1_pos, 2),
                              filename_neg = rep(scores_file_1_neg, 2),
                              filename_seg = rep(features_file_1, 2),
                              seqinfo = si)
  expect_type(sut, "list")
  expect_length(sut, 2)
  sut <- sut[[1]]
  expect_s4_class(sut, "NETseqData")
  total_reads <- sum(scores(sut)$score)
  expected_reads <- sum(scores$score * width(scores))
  expect_equal(total_reads, expected_reads)
  #TODO Better validity tests
})

test_that("NETseqData from GRanges scores", {
  nsd <-NETseqData(scores = unlist(tdata_gr_scores), sampleId = "S01")
  reads_from_granges <- sum(sapply(tdata_scores, function(u) 
    sum(unlist(u, use.names = FALSE))))
  reads_from_NETseqData <- sum(scores(nsd)$score)
  expect_equal(reads_from_NETseqData, reads_from_granges)
})


test_that("NETseqData from BAM file",{
  
  
  gene_list <- import("/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/genelist.gff", genome = "sacCer3")
  n_genes <- 20
  
  names(gene_list) <- gene_list$ID 
  
  z <- disjoin(gene_list, with.revmap = TRUE, ignore.strand = FALSE)$revmap
  w <- unique(unlist(z[which(sapply(z, function(u) length(u) > 1))]))
  
  if (length(w) > 0) {
    gene_list <- gene_list[-w]
  }
  
  if (n_genes > 0 && n_genes < length(gene_list)) {
    gene_list <- sample(gene_list, n_genes)
  }
  
  gene_list <- GenomicRanges::sort(gene_list)
  
  
  bam_directory <- "/n/groups/churchman/rds19/data/S005/mm-to-censor/"
  f <- tibble::tibble(sample_id = paste0("SRR1284006", 6:9),
              bam_file = paste0(bam_directory, sample_id, ".bam"))
  f <- f[1:2,]
  si <- Seqinfo(genome = "sacCer3")
  sut <- NETseqDataFromBAM(f$sample_id, f$bam_file, gene_list, si)
  
  expect_type(sut, "list")
  expect_true(all(sapply(sut, isa,"NETseqData")))
  # TODO: this should be in the sut
  names(sut) <- f$sample_id
  # TODO Break into 2 tests...after setting up appropriate fixture
  
  nsd <- sut
  sut <- GPosExperiment(rowRanges = gene_list, sample = nsd, seqinfo = si)
  
  expect_s4_class(sut, "GPosExperiment")
  # TODO: more verifcations
})
