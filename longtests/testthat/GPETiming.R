# GPETiming.R GPosExperiment timing
# TODO hardwired "real file names
f1 <- "/n/groups/churchman/rds19/data/S005/wt-1.pos.bedgraph.gz"
f2 <- "/n/groups/churchman/rds19/data/S005/wt-1.neg.bedgraph.gz"
frr <- "/n/groups/churchman/rds19/data/S005/genelist.gff"
rr <- import(frr)

si <- Seqinfo(genome = "sacCer3")

system.time(nsd <- NETseqDataFromBedgraph(sampleId = "wt-1", filename_pos = f1, filename_neg = f2, seqinfo = si))
system.time(sut <- GPosExperiment(sample = nsd, seqinfo = si, rowRanges = rr))
system.time(x <- scores(sut))
system.time(x <- scores(sut[1:10,]))
