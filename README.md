# GPosExperiment

### CAUTION: EXPERIMENTAL AND INCOMPLETE

GPosExperiment is an R package that extends Bioconductor's 
[SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
package.
The name GPosExperiment is derived from the GenomicRanges' GPos class, which is
designed to be especially efficient when representing values at the single
nucleotide resolution.
Its purpose to to facilitate the use of single-nucleotide resolution measurements 
in RangedSummarizedExperiment classes.

While motivated by the desire to exploit SummarizedExperiement containers to
manage [NET-seq](https://pubmed.ncbi.nlm.nih.gov/22470065/) datasets, it is 
intended to be useful for other assays which report values at
 genomic resolution.


