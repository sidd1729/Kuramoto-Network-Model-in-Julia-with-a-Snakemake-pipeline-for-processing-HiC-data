library(factoextra)
library(arrow)
library(cluster)
library(maotai)
library(dplyr)
library(pheatmap)
df <- read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/dcHiC_Regen_filtered_non_clustered_compartments.feather") #("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC1_100000.feather")#
#df_bin <- df$bin #When using dcHiC Output
df[df==-1]<-NA
df <- filter(df,!!sym(names(df)[1]) != -1 & !!sym(names(df)[2]) != -1 & !!sym(names(df)[3]) != -1 & !!sym(names(df)[4]) != -1 & !!sym(names(df)[5]) != -1) #When using cooltools output
df_omit <- na.omit(df)
df_choose <- names(df)[c(4,5,6,7,8)]
df <- df_omit[df_choose]
df_bin <- df$Chomosome_Region #When using cooltools output
#df_mat <- df[,c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")] #When using dcHiC output
#df_mat <- df[,c("control","12hpi","2dpi","4dpi","7dpi")]
df_mat <- as.matrix(df)
row.names(df_mat) <- df_bin
df_mat <- t(scale(t(df_mat)))


#create plot of number of clusters vs total within sum of squares
fviz_nbclust(df_mat, kmeans, method = "silhouette",k.max= 15)#"silhouette")#"wss")#,k.max=10) #"gap_stat"

starting_cluster_center <- kmeanspp(df_mat,k = 7)
empty_list = list(1:5)
full_list = rep(empty_list,7)
#empty_vec[[1]]
for (i in seq(1,7,by = 1)){
obj <- df[as.integer(names(starting_cluster_center[starting_cluster_center == i])),]
mean_obj <- obj %>% summarise(control_mean = mean(control_100kb),hpi_12_mean = mean(hpi_12_100kb), dpi_2_mean = mean(dpi_2_100kb),dpi_4_mean = mean(dpi_4_100kb),dpi_7_mean = mean(dpi_7_100kb))
full_list[[i]] = unlist(mean_obj[1,],use.names = FALSE) #c(mean_obj$control_mean,mean_obj$hpi_12_mean,mean_obj$dpi_2_mean,mean_obj$dpi_2_mean,mean_obj$dpi_4_mean,mean_obj$dpi_7_mean)
print(obj)
}
Check <- matrix(unlist(full_list),nrow = 7,ncol = 5)
k_mean_out<-kmeans(df_mat,centers = 7, nstart = 10000)
k_mean_out_center <- k_mean_out$centers
print(k_mean_out$centers)

mean <- colMeans(k_mean_out_center)
print(mean)

pheatmap(k_mean_out_center,cluster_row = FALSE, cluster_cols = FALSE)