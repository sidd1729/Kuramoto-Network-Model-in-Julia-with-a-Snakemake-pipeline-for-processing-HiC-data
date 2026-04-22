library(arrow)
library(pheatmap)
library(stringr)
library(dplyr)
library(furrr)
library(parallel)
library(tidyr)

nature <- "raw"#"raw"
sample_wise_df <- read_feather("/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.pcOri.feather")#"/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.Filtered.pcOri.feather")
sample_wise_df <- sample_wise_df %>% mutate(bin = str_glue("{chr}:{start}-{end}"))
rownames(sample_wise_df) <- sample_wise_df$bin
print(head(sample_wise_df))
#help(apply)
png(str_glue("results/dcHiC_Compartment_Cluster/dcHiC_Regen_filtered_compartments_{nature}.png"))#,res = 300)
len_unique_row <- apply(sample_wise_df,1,function (row) length(unique(row))) #Lambda function equivalent first taking unique elements of row and length of the row
#unique_row <- len_unique_row[len_unique_row != 1] #Identifying rows with unique elements
#print(dim(sample_wise_df))
#print(dim(sample_wise_df[unique_row,]))
df_to_hmap <- sample_wise_df[len_unique_row != 1,] #Choosing rows where len_unique_row !=1
df_hmap_bins <- df_to_hmap$bin
df_hmap_chr <- df_to_hmap$chr
df_hmap_start <- df_to_hmap$start
df_hmap_end <- df_to_hmap$end

df_to_hmap_2 <- df_to_hmap
#head(sample_wise_df)
df_to_hmap_2$bin_no <- rownames(sample_wise_df[len_unique_row != 1,])
write_feather(df_to_hmap_2,"results/dcHiC_Compartment_Cluster/dcHiC_Regen_filtered_non_clustered_compartments.feather")

df_to_hmap <- df_to_hmap[,c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")]
print(head(df_to_hmap))
#print(rownames(sample_wise_df[len_unique_row != 1,]))

#print(rownames(df_to_hmap))

#print(all(rownames(df_to_hmap) == rownames(sample_wise_df[len_unique_row != 1,])))
#rownames(df_to_hmap) <- rownames(sample_wise_df[len_unique_row != 1,])
#rownames(df_to_hmap) <- rownames(sample_wise_df[len_unique_row != 1,])
##df_to_hmap_2 <- df_to_hmap
#head(sample_wise_df)
##df_to_hmap_2$bin_no <- rownames(sample_wise_df[len_unique_row != 1,])
##write_feather(df_to_hmap_2,"results/dcHiC_Compartment_Cluster/dcHiC_Regen_filtered_non_clustered_compartments.feather")
#print(dim(df_to_hmap))
#x11()

# Quantile normalization portion
 if(nature == "quantiled"){
df_to_mat  <- as.matrix(df_to_hmap)
Q <-apply(apply(df_to_mat,2,rank), 2, function(r) mu=rowMeans(apply(df_to_mat, 2, sort))[r])
mat_dimensions <- dim(Q)
df_to_hmap[1:mat_dimensions[1],1:mat_dimensions[2]] <- Q}

#print(mat_dimensions)
res <- pheatmap(df_to_hmap,scale = "row",clustering_distance_rows="euclidean",cluster_cols = FALSE,show_rownames = FALSE)
res
invisible(dev.off())


sample_df_clust <- cbind(df_to_hmap, cluster = cutree(res$tree_row, k = 8))
sample_df_clust$rownames <- df_hmap_bins
sample_df_clust <-sample_df_clust[order(sample_df_clust$cluster),]
write_feather(sample_df_clust,str_glue("results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_{nature}_df.feather"))




print(head(sample_df_clust))

png(str_glue("results/dcHiC_Compartment_Cluster/dcHiC_Regen_filtered_compartments_{nature}_k_means.png"))#,res = 300)
#set.seed(123) Gives suboptimal results


df_mat <- df_to_hmap[,c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")]
row.names(df_mat) <- df_hmap_bins
#df_mat <- scale(df_mat) #Column scaling
df_mat <- t(scale(t(df_mat))) #Row scaling

k_mean_out<-kmeans(df_mat,centers = 8, nstart = 10000)
k_mean_out_center <- k_mean_out$centers
print(k_mean_out$centers)

mean <- colMeans(k_mean_out_center)
print(mean)

res <- k_mean_out
pheatmap_res <- pheatmap(k_mean_out_center,cluster_row = FALSE, cluster_cols = FALSE)
#res <- pheatmap(df_to_hmap,kmeans_k = 7, scale = "row",clustering_distance_rows="euclidean",cluster_cols = FALSE,show_rownames = TRUE)
invisible(dev.off())

sample_df_k_mean_clust = cbind(df_to_hmap,cluster = res$cluster)
sample_df_k_mean_clust$rownames <- df_hmap_bins
sample_df_k_mean_clust$chr <- df_hmap_chr
sample_df_k_mean_clust$start <- df_hmap_start
sample_df_k_mean_clust$end <- df_hmap_end
sample_df_k_mean_clust <- sample_df_k_mean_clust[order(sample_df_k_mean_clust$cluster),]

sample_df_k_mean_center <- data.frame(res$centers)

#print(sample_df_k_mean_center)
write_feather(sample_df_k_mean_clust,str_glue("results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_{nature}_k_means_df.feather"))
write_feather(sample_df_k_mean_center,str_glue("results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_{nature}_transformed_k_means_centers_df.feather"))
#print(sample_df_clust[order(sample_df_clust)])
#sSys.sleep(100)
print(head(sample_df_k_mean_clust,5))
plan(multisession,workers = 8)
cluster_centre_calc <- function(center_no,clustered_df = sample_df_k_mean_clust){#,column_names = c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")){
    print(names(clustered_df))
    #print(clustered_df[clustered_df$cluster == 1,])
    cluster_specific_df <- clustered_df[clustered_df$cluster == center_no,]
    #print(sym(column_names[1]))
    #cluster_specific_df <- cluster_specific_df %>% drop_na()
    #test_name <- as.name(column_names[1])
    #print(cluster_specific_df$test_name)
    mean_obj <- cluster_specific_df %>% summarise(across(c(control_100kb,hpi_12_100kb,dpi_2_100kb,dpi_4_100kb,dpi_7_100kb), ~mean(.x,na.rm = TRUE)))
    print(mean_obj)
    return (unlist(mean_obj[1,],use.names = FALSE))
}

mapped_example  <-future_map(seq(1,8,by = 1),cluster_centre_calc)
#print(mapped_example)
mat_test <-do.call("rbind",mapped_example)
#mat_test<-as.matrix(do.call("rbind",mapped_example))
#print(is.matrix(mat_test))
#print(is.matrix(mat_test_2))

k_means_unscaled <- data.frame(k_mean_out_center)
k_means_unscaled[1:dim(mat_test)[1],1:dim(mat_test)[2]] <- mat_test
png(str_glue("results/dcHiC_Compartment_Cluster/dcHiC_Regen_filtered_compartments_{nature}_unscaled_withclusteringscaled_k_means.png"))
pheatmap_res <- pheatmap(k_means_unscaled,cluster_row = FALSE, cluster_cols = FALSE)
#res <- pheatmap(df_to_hmap,kmeans_k = 7, scale = "row",clustering_distance_rows="euclidean",cluster_cols = FALSE,show_rownames = TRUE)
invisible(dev.off())
write_feather(k_means_unscaled,str_glue("results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_{nature}_k_means_centers_df.feather"))