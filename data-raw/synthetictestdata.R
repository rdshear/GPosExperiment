# Synthetic test data
# create short reads from selected sacCer3 genes
library(polyester)
library(Biostrings)
library(GenomicRanges)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# data_dir <- tempdir()
data_dir <- "~/Downloads"
chr <- c("chrI", "chrII")
g <- genes(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
Xcoverage <- 10
readlen <- 43
seed <- 20221120

setwd(data_dir)

si <- seqinfo(g)
r <- GRanges(seqname = chr, 
                 ranges =  rep(IRanges(start = 50000, end = 70000), length(chr)), 
                 seqinfo = seqinfo(g))

g <- g[subjectHits(findOverlaps(r,
                                g,
                                ignore.strand = TRUE))]

# TODO debug only
# shrink to 4 ROWs 
x <- length(g)
g <- g[c(1, 2, x-1, x)]
fa_filename <- "sacCer3_tx_fragment.fa"

seqs <- getSeq(BSgenome.Scerevisiae.UCSC.sacCer3,  GRanges(g, strand="*"))
writeXStringSet(seqs, fa_filename)
readspertx = round(Xcoverage * width(seqs) / 100)


fold_changes <- matrix(rep(1:2, times = length(seqs)), 
                       nrow = length(seqs), byrow = TRUE)

# TODO Distribution of read lengths?
# TODO: Use GTF + DNA sequence
# simulate_experiment(fa_filename, reads_per_transcript=readspertx, 
#                     paired=FALSE,
#                     seed=seed,
#                     readlen=readlen,
#                     num_reps=c(10,10), 
#                     fold_changes=fold_changes, 
#                     outdir='simulated_reads')

simulate_experiment(fa_filename, 
                    reads_per_transcript=readspertx, 
                    paired=FALSE,
                    seed=seed,
                    readlen=readlen,
                    num_reps=c(10,10),
                    fold_changes=fold_changes, 
                    outdir='simulated_reads')

