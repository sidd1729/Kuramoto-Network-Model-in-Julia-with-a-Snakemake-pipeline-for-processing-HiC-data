library(TADCompare)
library(stringr)
library(dplyr)
dataset <- "control_S68"
dataset_2 <- "12_hpi_S207"
dataset_3 <- "2_dpi_S208"
dataset_4 <- "4_dpi_S69"
dataset_5 <- "7_dpi_S1"

comparison_chrom <- 2
comparison_tp_1 <- 1
comparison_tp_2 <- 2
mats_tp_list <- c("sparse_mats","sparse_mats_2","sparse_mats_3","sparse_mats_4","sparse_mats_5")
tads_tp_list <- c("spec_tads","spec_tads_2","spec_tads_3","spec_tads_4","spec_tads_5")
#tested <- get(tads_tp_list[1])
#tested$chr1 == tested[[1]]
#tested[1]
sparse_mats = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset}_50kb_sparse.rds"))
sparse_mats_2 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_2}_50kb_sparse.rds"))
sparse_mats_3 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_3}_50kb_sparse.rds"))
sparse_mats_4 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_4}_50kb_sparse.rds"))
sparse_mats_5 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_5}_50kb_sparse.rds"))

spec_tads = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset}_50KB_TADs.rds"))
spec_tads_2 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_2}_50KB_TADs.rds"))
spec_tads_3 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_3}_50KB_TADs.rds"))
spec_tads_4 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_4}_50KB_TADs.rds"))
spec_tads_5 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_5}_50KB_TADs.rds"))

tested <- get(tads_tp_list[comparison_tp_1])
tested$chr1 == tested[[1]]
chr_bound_1 = bind_rows(get(tads_tp_list[comparison_tp_1])[[comparison_chrom]])
chr_bound_2 = bind_rows(get(tads_tp_list[comparison_tp_2])[[comparison_chrom]])
bound_bed <- list(chr_bound_1,chr_bound_2)
results = TADCompare(get(mats_tp_list[comparison_tp_1])[[comparison_chrom]],get(mats_tp_list[comparison_tp_2])[[comparison_chrom]],resolution = 50000, pre_tads = bound_bed)

###
#all(get(mats_tp_list[comparison_tp_1])$chr1 == get(mats_tp_list[comparison_tp_1])[[1]])
#all(get(mats_tp_list[comparison_tp_1])[[as.symbol(str_glue("chr{comparison_chrom}"))]] == sparse_mats$chr1)
###

# Visualizing the results
p <- DiffPlot(tad_diff  = results, 
              cont_mat1   = get(mats_tp_list[comparison_tp_1])[[comparison_chrom]],
              cont_mat2   = get(mats_tp_list[comparison_tp_2])[[comparison_chrom]],
              resolution  = 50000,
              start_coord = 1,
              end_coord   = 6500000,
              pre_tad     = bound_bed,
              show_types  = FALSE, 
              point_size  = 3,
              palette     = "Spectral",
              rel_heights = c(1, 1))
plot(p)