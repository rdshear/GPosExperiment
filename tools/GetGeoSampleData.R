# Get sample data for HresExperiment

library(GEOquery)
library(BSgenome.Hsapiens.UCSC.hg38)

# GSE223034 GSE223036

GSE223034 <- getGEOSuppFiles("GSE223034", baseDir = "~/temp", fetch_files = TRUE)
x <- getGEOfile("GSE223034", amount = "full")
y <- getGSEDataTables("GSE223034")
z <- getGEO("PRJNA924743", destdir = "~/temp")
