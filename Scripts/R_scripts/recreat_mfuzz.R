#Recreating mfuzz plot from centers and membership matrix
library(Mfuzz)
library(arrow)
library(stringr)
library(furrr)
library(dplyr)
library(pheatmap)
library(factoextra)
nature <- "raw"
df_type <- "dcHiC_PC"
image_type <- "svg"

#To be commented out later
sample_wise_df <- read_feather("/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.pcOri.feather")#"/home/ksslab/Programs/dcHiC/DifferentialResult/Regen_100kb_timepoint/fdr_result/differential.intra_sample_group.Filtered.pcOri.feather")
sample_wise_df <- sample_wise_df %>% mutate(bin = str_glue("{chr}:{start}-{end}"))
rownames(sample_wise_df) <- sample_wise_df$bin
##print(head(sample_wise_df))
#help(apply)
len_unique_row <- apply(sample_wise_df,1,function (row) length(unique(row))) #Lambda function equivalent first taking unique elements of row and length of the row
df_to_hmap <- sample_wise_df[len_unique_row != 1,] #Choosing rows where len_unique_row !=1
df_hmap_bins <- df_to_hmap$bin
df_hmap_chr <- df_to_hmap$chr
df_hmap_start <- df_to_hmap$start
df_hmap_end <- df_to_hmap$end

df_to_hmap_2 <- df_to_hmap
###head(sample_wise_df)
df_to_hmap_2$bin_no <- rownames(sample_wise_df[len_unique_row != 1,])
#write_feather(df_to_hmap_2,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_Regen_filtered_non_clustered_compartments_yan.feather"))

##df_to_hmap <- df_to_hmap[,c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")]
DLR_df <- df_to_hmap
DLR_choose = names(DLR_df[c(4,5,6,7,8)])
###

DLR_df_omit <- na.omit(DLR_df)
DLR_filtered <- DLR_df_omit[DLR_choose]
M <- as.matrix(DLR_filtered)
if(nature == "quantiled"){
  Q=apply(apply(M,2,rank), 2, function(r) mu=rowMeans(apply(M, 2, sort))[r])} else{
    Q <- M}
eset <- new('ExpressionSet', exprs=Q)
standardized <- standardise(eset)
m1 <- mestimate(standardized)

DLR_df_bound <-read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_clustered_cooltools_PCs_sample_wise_{nature}_k_means_df_yan.feather"))
plan(multisession,workers = 8)
cluster_centre_calc <- function(center_no,clustered_df = DLR_df_bound){#,column_names = c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")){
  print(names(clustered_df))
  #print(clustered_df[clustered_df$cluster == 1,])
  cluster_specific_df <- clustered_df[clustered_df$cluster == center_no,]
  #print(sym(column_names[1]))
  #cluster_specific_df <- cluster_specific_df %>% drop_na()
  #test_name <- as.name(column_names[1])
  #print(cluster_specific_df$test_name)
  #1st line below for cooltools, 2nd on for dchic,3rd for K_50
  #mean_obj <- cluster_specific_df %>% summarise(across(c(control,`12hpi`,`2dpi`,`4dpi`,`7dpi`), ~mean(.x,na.rm = TRUE)))
  mean_obj <- cluster_specific_df %>% summarise(across(c(`control_100kb`,`hpi_12_100kb`,`dpi_2_100kb`,`dpi_4_100kb`,`dpi_7_100kb`), ~mean(.x,na.rm = TRUE)))
  #mean_obj <- cluster_specific_df %>% summarise(across(c(control,`12 hpi`,`2 dpi`,`4 dpi`,`7 dpi`), ~mean(.x,na.rm = TRUE)))
  print(mean_obj)
  return (unlist(mean_obj[1,],use.names = FALSE))
}

mapped_example  <-future_map(seq(1,8,by = 1),cluster_centre_calc)
#print(mapped_example)

mat_test <-do.call("rbind",mapped_example)
names(mat_test)
df_test <- as.data.frame(mat_test)
names(df_test) <- c("control","12hpi","2dpi","4dpi","7dpi")
##png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_Cooltools_Regen_compartments_{nature}_unscaled_withclusteringscaled_k_means_yan.png"))
svg(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_Cooltools_Regen_compartments_{nature}_unscaled_withclusteringscaled_k_means_yan.svg"))
pheatmap_res <- pheatmap(df_test,cluster_row = FALSE, cluster_cols = FALSE)
pheatmap_res
dev.off()

#-------
# Create the object
centroid_df <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_clustered_cooltools_PCs_sample_wise_{nature}_k_means_centers_df_yan.feather"))
membership_df <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_clustered_cooltools_PCs_sample_wise_{nature}_transformed_k_means_membership_df_yan.feather"))
cluster <- DLR_df_bound$cluster
mem_mat <- as.matrix(membership_df)
cluster_vector <- apply(mem_mat, 1, which.max)


cl_object <- list(
  as.matrix(centroid_df),        
  NA,              
  cluster_vector,  
  mem_mat          
)
if (image_type == "png"){
  png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_{nature}_yan.png"))
} else if (image_type == "svg"){
  svg(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_{nature}_yan.svg"))}
par(mar = c(5, 4, 4, 1) + 0.1, oma = c(0, 0, 3, 0))
mfuzz.plot2(standardized,cl=cl_object,mfrow=c(4,2),time.labels=c("control","12hpi","2dpi","4dpi","7dpi"),x11 = FALSE,centre = FALSE)#,#centre.col = "blue")

title(str_glue("{df_type} clusters"), outer = TRUE, cex.main = 2.0)
invisible(dev.off())

my_pca <- prcomp(centroid_df, scale. = TRUE, center = TRUE, retx = TRUE)
#names(my_pca)
if (image_type == "png"){
  png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/PCA/Contributions/Contributions_{df_type}_{nature}_yan.png"))
} else if (image_type == "svg"){
  svg(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/PCA/Contributions/Contributions_{df_type}_{nature}_yan.svg"))}
fviz_pca_var(my_pca,title = str_glue("{df_type} {nature} contributions"))
invisible(dev.off())

if (image_type == "png"){
  png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/PCA/Cluster_centroids/Cluster_centroids_{df_type}_{nature}_yan.png"))
} else if (image_type == "svg"){
  svg(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/PCA/Cluster_centroids/Cluster_centroids_{df_type}_{nature}_yan.svg"))}
fviz_pca_ind(my_pca, title = str_glue("{df_type} {nature} cluster centroids"))
invisible(dev.off())