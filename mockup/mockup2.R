library(rtracklayer)

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
}

library(BSgenome.Hsapiens.UCSC.hg38)
gn <- BSgenome.Hsapiens.UCSC.hg38

setwd("~/temp")
x <- rtracklayer::import.wig("GSM6938441_caps_4iesc_1.bw")
head(x)
y <- import.wig("GSM6938447_CR_d8DZPGCLC_DMRT1_r1.bigWig", genome = 'hg38')

a <- GRangesList(lapply(list(x,y), (\(u) u[seqnames(u) == "chr1"])))
ov <- findOverlaps(a[[1]], a[[2]])

lapply(a, (\(u) length(ov) / length(u)))
