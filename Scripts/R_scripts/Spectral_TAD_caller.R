#Hierarchical (or not) TAD calling with spectral TAD
library(SpectralTAD)
library(stringr)
library(purrr)
#cool_mat = read.table("~/Downloads/Rao.GM12878.50kb.txt")


#data("rao_chr20_25_rep")
#head(rao_chr20_25_rep)

#We see that this is a sparse 3-column contact matrix
#Running the algorithm with resolution specified
#results = SpectralTAD(rao_chr20_25_rep, chr = "chr20", resolution = 25000, qual_filter = FALSE, z_clust = FALSE)
#Printing the top 5 TADs
#head(results$Level_1, 5)

#Repeating without specifying resolution
#no_res = SpectralTAD(rao_chr20_25_rep, chr = "chr20", qual_filter = FALSE, z_clust = FALSE)
#We can see below that resolution can be estimated automatically if necessary
#identical(results, no_res)

for (file in list("control_S68_50kb")){
args <- commandArgs(trailingOnly = TRUE)
cool_mat = read.table(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_txt/{file}.txt"))
cool_mat$V1 <- unlist(map(cool_mat$V1,function(x) str_glue("chr{x}")))
cool_mat$V4 <- unlist(map(cool_mat$V4,function(x) str_glue("chr{x}")))
#Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
sparse_mats = HiCcompare::cooler2sparse(cool_mat)
    #spec_tads = lapply(names(sparse_mats), function(x) {
  #SpectralTAD(sparse_mats[[x]], chr = x)
#})



saveRDS(sparse_mats,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{file}_sparse.rds"))
readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{file}_sparse.rds"))}