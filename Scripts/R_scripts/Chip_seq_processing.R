#Processing alignment files, bam to give reads in a given bin hopefully?
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(Rsamtools)
filename = "C:/Users/HP/Downloads/ENCFF582IQY.bam"
bf <- BamFile(filename)
bf
seqinfo(bf)
scanBam(BamFile(filename, yieldSize=5))
import.bed("C:/Users/HP/Downloads/ENCFF582IQY.bam")

which <- GRanges(c(
  "chr1:1000-2000"
))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which)#, what=what)
bamFile <- BamFile(filename)
bam <- scanBam(bamFile, param=param)
header <- scanBamHeader(bf)
filename = "~/Siddharth/Program_Directory/Data/ChiP_Seq_Data/H1_HESC/CTCF/ENCFF616PLK.bam"
Index_bam <- indexBam(filename)
bf = open(BamFile(filename, Index_bam))
bf