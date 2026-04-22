library(arrow)
mem_mat <- read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_raw_transformed_k_means_membership_df_yan.feather")

#Plotting frequency histograms if needed
#hist(mem_mat$X8,xlab ="Cluster 8 membership",ylab = "Frequency of Cluster 8",main = "Distribution of cluster 8")
cluster_1_den <- density(mem_mat$X1)
cluster_2_den <- density(mem_mat$X2)
cluster_3_den <- density(mem_mat$X3)
cluster_4_den <- density(mem_mat$X4)
cluster_5_den <- density(mem_mat$X5)
cluster_6_den <- density(mem_mat$X6)
cluster_7_den <- density(mem_mat$X7)
cluster_8_den <- density(mem_mat$X8)
overall_dens <- density(c(mem_mat$X1,mem_mat$X2,mem_mat$X3,mem_mat$X4,mem_mat$X5,mem_mat$X6,mem_mat$X7,mem_mat$X8))

#Making density plots
png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_1_membership_density.png")
plot(cluster_1_den$x,cluster_1_den$y,type = "l",xlab ="Cluster 1 membership",ylab = "Density of Cluster 1")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_2_membership_density.png")
plot(cluster_2_den$x,cluster_2_den$y,type = "l",xlab ="Cluster 2 membership",ylab = "Density of Cluster 2")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_3_membership_density.png")
plot(cluster_3_den$x,cluster_3_den$y,type = "l",xlab ="Cluster 3 membership",ylab = "Density of Cluster 3")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_4_membership_density.png")
plot(cluster_4_den$x,cluster_4_den$y,type = "l",xlab ="Cluster 4 membership",ylab = "Density of Cluster 4")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_5_membership_density.png")
plot(cluster_5_den$x,cluster_5_den$y,type = "l",xlab ="Cluster 5 membership",ylab = "Density of Cluster 5")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_6_membership_density.png")
plot(cluster_6_den$x,cluster_6_den$y,type = "l",xlab ="Cluster 6 membership",ylab = "Density of Cluster 6")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_7_membership_density.png")
plot(cluster_7_den$x,cluster_7_den$y,type = "l",xlab ="Cluster 7 membership",ylab = "Density of Cluster 7")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_cluster_8_membership_density.png")
plot(cluster_8_den$x,cluster_8_den$y,type = "l",xlab ="Cluster 8 membership",ylab = "Density of Cluster 8")#,cex.axis = 0.5)
dev.off()

png("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_all_clusters_membership_density.png")
plot(cluster_1_den$x,cluster_1_den$y,type = "l",xlab ="All clusters membership",ylab = "Density of all clusters")#,cex.axis = 0.5)
dev.off()

#Printing the general info for all
print(cluster_1_den)
print(cluster_2_den)
print(cluster_3_den)
print(cluster_4_den)
print(cluster_5_den)
print(cluster_6_den)
print(cluster_7_den)
print(cluster_8_den)
print(overall_dens)
