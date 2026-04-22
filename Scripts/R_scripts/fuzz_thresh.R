library(arrow)
library(stringr)
library(furrr)
library(pheatmap)
df_type <- "dcHiC_PC" #"K_50"
image_type <- "svg"
nature <- "raw"
clus_1_thres <- 0.33320
clus_2_thres <- 0.26252
clus_3_thres <- 0.537021
clus_4_thres <- 0.2592
clus_5_thres <- 0.3314
clus_6_thres <- 0.34474
clus_7_thres <- 0.33028
clus_8_thres <- 0.524386
clus_thres_vec <- c(clus_1_thres,clus_2_thres,clus_3_thres,clus_4_thres,clus_5_thres,clus_6_thres,clus_7_thres,clus_8_thres)
thres_decided = min(clus_thres_vec)

mem_mat <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_clustered_cooltools_PCs_sample_wise_{nature}_transformed_k_means_membership_df_yan.feather"))
write.csv(mem_mat,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/{df_type}_clustered_cooltools_PCs_sample_wise_{nature}_transformed_k_means_membership_df_yan.csv"))
#max_membership <- apply(mem_mat,1,max)
#max_40 <- max_membership[max_membership > 0.40]
print(mem_mat)
#print(max_membership)
#print(max_40)
thresh_assign <- function(row_vector,thres=thres_decided){
  appropriate_indices = which(row_vector >= thres)
  #print(appropriate_indices)
  #print(row_vector)
  #print(length(appropriate_indices))
  if (length(appropriate_indices)>1){
     membership_candidates = row_vector[appropriate_indices]
     appropriate_index = which(row_vector == max(membership_candidates))
     if(row_vector[appropriate_index] >= clus_thres_vec[appropriate_index]){
       return (appropriate_index)
    }
    else{
      return(0)
    }
  }
  else if (length(appropriate_indices) == 0){
    return(0)
  }
  else{
    if(row_vector[appropriate_indices] >= clus_thres_vec[appropriate_indices]){
     return (appropriate_indices)
  }
   else{
     return(0)
   }
  }
  }

row_helper <- function(i){
  return(thresh_assign(mem_mat[i,]))
}

#print(mem_mat[4,])
plan(multisession(workers = 8))
thresholded_cluster_labels <- unlist(future_map(seq(1,dim(mem_mat)[1],by = 1),row_helper))

dcHiC_df <- read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_Regen_filtered_non_clustered_compartments_yan.feather")
write.csv(dcHiC_df,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_Regen_filtered_non_clustered_compartments_yan.csv")
dcHiC_clust <- dcHiC_df %>% mutate(cluster = thresholded_cluster_labels)
write_feather(dcHiC_clust,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_compartments_yan.feather"))

cluster_centre_calc <- function(center_no,clustered_df = dcHiC_clust){#,column_names = c("control_100kb","hpi_12_100kb","dpi_2_100kb","dpi_4_100kb","dpi_7_100kb")){
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
png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_{nature}_thres_{thres_decided}_unscaled_withclusteringscaled_k_means_yan.png"))
pheatmap_res <- pheatmap(df_test,cluster_row = FALSE, cluster_cols = FALSE)
pheatmap_res
dev.off()
write_feather(df_test,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_{nature}_thres_{thres_decided}_k_means_centers_df_yan.feather"))