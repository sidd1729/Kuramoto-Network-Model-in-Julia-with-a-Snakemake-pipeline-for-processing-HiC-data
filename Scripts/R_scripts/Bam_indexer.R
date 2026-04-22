#Bam indice generator
args <- commandArgs(trailingOnly = TRUE)
library(Rsamtools)
indexBam(args[1])