# mock up prototype for assay approach
library(HRseqExperiment)

data("test_bedgraphs")
data("test_masks")
data("genelist")

genelist <-  genelist[1:5,c("curie", "gene")]
mcols(genelist) <- mcols$gene
sampleId = "SRR12840066"
refmasks <- test_masks$SRR12840066
nsd <- lapply(test_bedgraphs, function(u) {
  HRseqData(scores = u,
            mask = refmasks,
            sampleId = sampleId)
})
sut <- HRseqExperiment(sample = nsd, rowRanges = genelist,
                       seqinfo = refdata$seqinfo)

ref_scores <- .OverlappedRanges(genelist, test_bedgraphs$SRR12840066)
ref_totalreads <- sum(ref_scores$score * width(ref_scores))
s <- scores(sut, zero_fill = FALSE, apply_mask = FALSE)
sl <- apply(s, 1:2, function(u) length(u[[1]]))
sc <- apply(s, 1:2, function(u) sum(u[[1]]$score))
sz <- scores(sut, zero_fill = TRUE, apply_mask = FALSE)
szl <- apply(sz, 1:2, function(u) length(u[[1]]))
szc <- apply(sz, 1:2, function(u) sum(u[[1]]$score))
m <- mask(sut)
smat <- smat <-apply(sz, 1:2, function(u) u[[1]]$score)


grl <- test_bedgraphs
sample_names <- names(grl)

# HRseqAssay(x = GRangesList(), seqFilenameList = NULL, 
#            strand.aware = FALSE, zero.fill = FALSE, mcol.name = "score", value.type = numeric)
# 
# assay(x, i) <- HRseqAssay(...)
# 
# assay(se, i = 1, formula = ., return.type = numeric, strand.aware = TRUE)

genelist
str(sample_names)
str(grl, max.level = 1)
se <- SummarizedExperiment(rowRanges = genelist, colData = sample_names,
                           assays = list(score = smat))
# 
print(assay(se, score))
print(apply(s, 1:2, function(u) sum(u[[1]]$score)))
