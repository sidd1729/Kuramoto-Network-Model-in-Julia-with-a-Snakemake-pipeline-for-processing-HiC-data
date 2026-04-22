library(arrow)
library(stringr)
library(pheatmap)
#library(fclust)


df_to_hmap <- read_feather("results/k_50_table/filtered_non_clustered_k_50_sample_wise_df.feather")
row_abs_diff_sum <- apply(df_to_hmap,1,function (row) abs(row[1]-row[2]) + abs(row[2]-row[3]) + abs(row[3]-row[4]) + abs(row[4]-row[5]))
print(row_abs_diff_sum)

res <- pheatmap(row_abs_diff_sum,scale = "column",clustering_distance_rows="euclidean",cluster_rows = FALSE ,cluster_cols = FALSE ,show_rownames = FALSE)
res