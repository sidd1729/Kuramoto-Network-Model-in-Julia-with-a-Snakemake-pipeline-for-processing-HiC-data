library(arrow)
library(stringr)
library(fclust)

df_to_hmap <- read_feather("results/k_50_table/filtered_non_clustered_k_50_sample_wise_df.feather")
#fkm <- FKM(X = df_to_hmap, m = 1.2, RS = 50, stand = 1, index = "SIL.F")
#summary(fkm)

library(cluster)

gap_stat <- clusGap(df_to_hmap, FUN = kmeans, nstart = 25,
                    K.max = 25, B = 50)

print(gap_stat)

#fviz_gap_stat(gap_stat)