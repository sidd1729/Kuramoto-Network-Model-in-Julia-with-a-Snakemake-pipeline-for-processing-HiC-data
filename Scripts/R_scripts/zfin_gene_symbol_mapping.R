#Mapping zfin ID back to GO symbol if required
Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')
library(stringr)
library(enrichR)
library(dplyr)
library(ggplot2)
listEnrichrSites()

#Setting site as FishEnrichr
setEnrichrSite("FishEnrichr")

listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018","GO_Cellular_Component_2018","GO_Biological_Process_2018")
nature <- "raw"
thres_decided = 0.2592
zfin_map_req <- TRUE
zfin_annotation_file <- read.csv("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin.gaf",sep = '\t',comment.char = "!",header = FALSE)
names(zfin_annotation_file) <- c("Database_Designation","Marker_ID","Gene_Symbol","Qualifiers","GO_Term_ID",
                          "Reference_ID","GO_Evidence_Code","Inferred_From","Ontology","Marker_Name","Marker_Synonyms",
                          "Marker_Type","Taxon","Modification_Date","Assigned_By","Annotation_Extension","Gene_Product_Form_ID")

zfin_annotation_file_bp = zfin_annotation_file[zfin_annotation_file$Ontology == "P",]
#file_to_map <- read.csv("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_raw_k_means_cluster_all_genes.csv")
file_to_map <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_all_genes.csv"))
Required  <-inner_join(file_to_map,zfin_annotation_file_bp, by = join_by(x == Marker_ID))
#write.csv(unique(Required$Gene_Symbol),"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_raw_k_means_cluster_all_genes_symbols.csv")
write.csv(unique(Required$Gene_Symbol),str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_all_genes_symbols.csv"))

compartment_background_gene_set <- unique(Required$Gene_Symbol)

for (i in seq(1,8,by =1)){
#file_to_map <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{i}_genes.csv"))
file_to_map <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{i}_genes.csv"))
Required  <-inner_join(file_to_map,zfin_annotation_file_bp, by = join_by(x == Marker_ID))
#write.csv(unique(Required$Gene_Symbol),str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{i}_genes_symbols.csv"))
write.csv(unique(Required$Gene_Symbol),str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{i}_genes_symbols.csv"))
df <- enrichr(unique(Required$Gene_Symbol), databases = dbs, background = NULL, include_overlap = FALSE, sleepTime = 1)
#write.csv(df$GO_Biological_Process_2018,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_cluster_cool_Overall_regeneration_{i}_fishEnrichR_BP_Ontology.csv"))
write.csv(df$GO_Biological_Process_2018,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{i}_fishEnrichR_BP_Ontology.csv"))
#png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_cluster_cool_Overall_regeneration_{cat}_fishEnrichR_BP_Ontology.png"))

tempPlot1 <- plotEnrich(df$GO_Biological_Process_2018, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", xlab = NULL, ylab = NULL, title = str_glue("Cluster {i} BP"))
#ggsave(tempPlot1,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_cluster_cool_Overall_regeneration_{i}_fishEnrichR_BP_Ontology.png"))
ggsave(tempPlot1,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{i}_fishEnrichR_BP_Ontology.png"))
#dev.off()
}


